#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <ctime>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "Psapi.lib")
#else
#include <sys/resource.h>
#endif

using namespace std;

template<class K, class V>
using HashMap = unordered_map<K, V>;

template<class K>
using HashSet = unordered_set<K>;

// 염기 유효성 / 대문자 변환 / 인코딩(A=0,C=1,G=2,T=3) 룩업 테이블
static bool DNA_VALID[256];
static char DNA_UPPER[256];
static uint8_t BASE_ENC[256];
static const char BASE_DEC[4] = { 'A', 'C', 'G', 'T' };

// 룩업 테이블 초기화 — 프로그램 시작 시 1회 호출
static void init_tables()
{
	for (int i = 0; i < 256; i++) {
		DNA_VALID[i] = false;
		DNA_UPPER[i] = (char)i;
		BASE_ENC[i] = 0;
	}

	DNA_VALID['A'] = DNA_VALID['C'] = DNA_VALID['G'] = DNA_VALID['T'] = true;
	DNA_VALID['a'] = DNA_VALID['c'] = DNA_VALID['g'] = DNA_VALID['t'] = true;

	DNA_UPPER['a'] = 'A'; DNA_UPPER['c'] = 'C';
	DNA_UPPER['g'] = 'G'; DNA_UPPER['t'] = 'T';

	BASE_ENC['A'] = BASE_ENC['a'] = 0; BASE_ENC['C'] = BASE_ENC['c'] = 1;
	BASE_ENC['G'] = BASE_ENC['g'] = 2; BASE_ENC['T'] = BASE_ENC['t'] = 3;
}

// 현재 프로세스 메모리 사용량(MB) 반환 — OS별 분기
double get_memory_usage_mb()
{
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS pmc;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
		return (double)pmc.WorkingSetSize / (1024.0 * 1024.0);
	return 0.0;
#else
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage) == 0)
		return (double)usage.ru_maxrss / 1024.0;
	return 0.0;
#endif
}

// 염기 1개를 2비트 정수로 인코딩 / 디코딩
inline uint64_t encode_base(char c) { return BASE_ENC[(unsigned char)c]; }
inline char     decode_base(uint64_t x) { return BASE_DEC[x & 3]; }

// 2비트 인코딩된 k-mer 값을 문자열로 복원
string decode_kmer(uint64_t value, int length)
{
	string s(length, 'A');
	for (int i = length - 1; i >= 0; i--) { s[i] = decode_base(value); value >>= 2; }
	return s;
}

// FASTA 파일에서 염기서열만 추출 (헤더 '>' 줄 무시, max_len 지정 시 조기 종료)
string load_fasta(const string& filename, size_t max_len = 0)
{
	ifstream f(filename);
	if (!f.is_open()) { cerr << "FASTA open fail: " << filename << "\n"; return ""; }

	string genome, line;
	if (max_len > 0) genome.reserve(max_len);

	while (getline(f, line)) {
		if (line.empty() || line[0] == '>') continue;
		for (unsigned char c : line) {
			if (DNA_VALID[c]) {
				genome += DNA_UPPER[c];
				if (max_len > 0 && genome.size() >= max_len) {
					genome.shrink_to_fit(); return genome;
				}
			}
		}
	}
	genome.shrink_to_fit();
	return genome;
}

// 서열의 역상보(reverse complement) 반환
string reverse_complement(const string& s)
{
	string rc(s.size(), 'A');
	for (size_t i = 0; i < s.size(); i++) {
		char c = s[s.size() - 1 - i];
		if (c == 'A' || c == 'a') rc[i] = 'T';
		else if (c == 'T' || c == 't') rc[i] = 'A';
		else if (c == 'C' || c == 'c') rc[i] = 'G';
		else if (c == 'G' || c == 'g') rc[i] = 'C';
		else rc[i] = 'N';
	}
	return rc;
}

// 리드에 무작위 염기 치환(오류) 삽입 — 시퀀싱 오류 시뮬레이션
void add_random_mismatch(string& read, int max_mismatch)
{
	if (max_mismatch <= 0 || read.empty()) return;
	int num_mutations = rand() % (max_mismatch + 1);
	static const char bases[4] = { 'A', 'C', 'G', 'T' };
	for (int m = 0; m < num_mutations; m++) {
		int pos = rand() % (int)read.size();
		char original = read[pos], mutated;
		do { mutated = bases[rand() % 4]; } while (mutated == original);
		read[pos] = mutated;
	}
}

// paired-end 리드 한 쌍을 담는 구조체
// read2는 역상보 형태로 저장 (실제 시퀀서 출력 방향 반영)
struct PairedRead {
	string read1;
	string read2_rc;   // read2의 역상보
	int frag_start;
	int insert_size;
};

// 게놈에서 paired-end 리드 시뮬레이션
// 무작위 위치에서 fragment를 뽑아 read1/read2 생성 후 오류 삽입
vector<PairedRead> generate_paired_reads(
	const string& genome, int pair_count, int read_len,
	int insert_min, int insert_max, int max_mismatch)
{
	vector<PairedRead> pairs;
	pairs.reserve(pair_count);
	int n = (int)genome.size();

	for (int i = 0; i < pair_count; i++) {
		int frag_len = insert_min + rand() % (insert_max - insert_min);
		if (frag_len + read_len >= n) continue;

		int frag_start = rand() % (n - frag_len - read_len);

		string read1 = genome.substr(frag_start, read_len);
		add_random_mismatch(read1, max_mismatch);

		// read2는 fragment 끝에서 역방향으로 읽으므로 역상보 저장
		string read2 = genome.substr(frag_start + frag_len, read_len);
		string read2_rc = reverse_complement(read2);
		add_random_mismatch(read2_rc, max_mismatch);

		PairedRead pr;
		pr.read1 = read1; pr.read2_rc = read2_rc;
		pr.frag_start = frag_start;
		pr.insert_size = frag_len + read_len;
		pairs.push_back(pr);
	}
	return pairs;
}

// De Bruijn 그래프 간선 (목적지 노드 + 가중치)
struct Edge { uint64_t to; uint32_t weight; };

// De Bruijn 그래프 — 노드는 (K-1)-mer, 간선은 K-mer
class DeBruijnGraph
{
public:
	HashMap<uint64_t, vector<Edge>> adj;  // 인접 리스트
	int K;
	uint64_t node_mask;  // (K-1)-mer 비트 마스크

	explicit DeBruijnGraph(int k) : K(k)
	{
		assert(k >= 2 && k <= 32);
		node_mask = (K - 1 < 32) ? (1ULL << (2 * (K - 1))) - 1ULL : ~0ULL;
		adj.reserve(1 << 20);
	}

	// 리드 한 개를 그래프에 추가 — 슬라이딩 윈도우로 K-mer 간선 생성
	void add_read(const string& read)
	{
		int n = (int)read.size();
		if (n < K) return;

		// 첫 (K-1)-mer를 prefix로 초기화
		uint64_t prefix = 0;
		for (int i = 0; i < K - 1; i++)
			prefix = ((prefix << 2) & node_mask) | encode_base(read[i]);
		if (adj.find(prefix) == adj.end()) adj[prefix] = vector<Edge>();

		// 나머지 염기를 순회하며 prefix→suffix 간선 추가(중복이면 weight++)
		for (int i = K - 1; i < n; i++) {
			uint64_t b = encode_base(read[i]);
			uint64_t suffix = ((prefix << 2) & node_mask) | b;

			vector<Edge>& edges = adj[prefix];
			bool found = false;
			for (size_t j = 0; j < edges.size(); j++) {
				if (edges[j].to == suffix) { edges[j].weight++; found = true; break; }
			}
			if (!found) {
				Edge e; e.to = suffix; e.weight = 1;
				edges.push_back(e);
				if (adj.find(suffix) == adj.end()) adj[suffix] = vector<Edge>();
			}
			prefix = suffix;
		}
	}

	// 낮은 가중치(= 시퀀싱 오류성 간선) 제거 — 반환값: 제거된 간선 수
	int prune(int min_weight)
	{
		int removed_edges = 0;
		for (auto& kv : adj) {
			vector<Edge>& edges = kv.second;
			size_t before = edges.size();
			edges.erase(
				remove_if(edges.begin(), edges.end(),
					[min_weight](const Edge& e) { return (int)e.weight < min_weight; }),
				edges.end());
			removed_edges += (int)(before - edges.size());
		}
		return removed_edges;
	}
};

// ── [추가] paired-end guided assembly용 링크 테이블 ──────────────────
// PELinkTable[kmer_A][kmer_B] = 두 k-mer가 같은 리드 쌍에 속하는 횟수
using PELinkTable = HashMap<uint64_t, HashMap<uint64_t, int>>;

// 리드 쌍에서 k-mer 레벨 paired-end 링크 테이블 구축
// read1의 마지막 k-mer ↔ read2의 첫 k-mer 를 양방향으로 기록
PELinkTable build_pe_link_table(const vector<PairedRead>& pairs, int k)
{
	PELinkTable table;
	table.reserve(1 << 20);
	uint64_t mask = (k < 32) ? (1ULL << (2 * k)) - 1ULL : ~0ULL;

	auto extract_kmers = [&](const string& read) -> vector<uint64_t> {
		vector<uint64_t> kmers;
		if ((int)read.size() < k) return kmers;
		uint64_t val = 0;
		for (int i = 0; i < k - 1; i++)
			val = ((val << 2) & mask) | encode_base(read[i]);
		for (int i = k - 1; i < (int)read.size(); i++) {
			val = ((val << 2) & mask) | encode_base(read[i]);
			kmers.push_back(val);
		}
		return kmers;
	};

	for (size_t i = 0; i < pairs.size(); i++) {
		string r2 = reverse_complement(pairs[i].read2_rc);
		vector<uint64_t> kmers1 = extract_kmers(pairs[i].read1);
		vector<uint64_t> kmers2 = extract_kmers(r2);

		// read1 끝 N개 k-mer ↔ read2 앞 N개 k-mer
		const int ANCHOR = 5;

		vector<uint64_t> k1_anchors(
			kmers1.size() >= ANCHOR ? kmers1.end() - ANCHOR : kmers1.begin(),
			kmers1.end());
		vector<uint64_t> k2_anchors(
			kmers2.begin(),
			kmers2.size() >= ANCHOR ? kmers2.begin() + ANCHOR : kmers2.end());

		for (uint64_t k1 : k1_anchors) {
			for (uint64_t k2 : k2_anchors) {
				table[k1][k2]++;
				table[k2][k1]++;
			}
		}
	}
	return table;
}

// unitig 방식으로 contig 추출
// 분기 없는 직선 경로(in=1, out=1)를 이어 붙여 하나의 contig로 만듦
// 분기점에서 paired-end 링크 테이블로 경로 선택
vector<string> generate_contigs_unitig(DeBruijnGraph& dbg, int min_len,
	const PELinkTable& pe_links)
{
	// 각 노드의 입력/출력 차수 계산
	HashMap<uint64_t, int> indeg, outdeg;
	for (auto& kv : dbg.adj) {
		uint64_t from = kv.first;
		outdeg[from] = (int)kv.second.size();
		if (indeg.find(from) == indeg.end()) indeg[from] = 0;
		for (size_t i = 0; i < kv.second.size(); i++) indeg[kv.second[i].to]++;
	}

	// 간선 방문 여부를 유일 키로 추적
	auto edge_key = [](uint64_t from, uint64_t to) -> uint64_t {
		return from ^ (to * 2654435761ULL);
	};

	// 분기점 또는 시작점(in≠1 또는 out≠1) 여부 판단
	auto is_start = [&](uint64_t v) -> bool {
		int id = indeg.count(v) ? indeg[v] : 0;
		int od = outdeg.count(v) ? outdeg[v] : 0;
		if (od == 0) return false;
		return !(id == 1 && od == 1);
	};

	HashSet<uint64_t> visited;
	vector<string> contigs;

	// pass 0: 분기점에서 출발 / pass 1: 미방문 간선 처리
	for (int pass = 0; pass < 2; pass++) {
		for (auto& kv : dbg.adj) {
			uint64_t start = kv.first;
			if (pass == 0 && !is_start(start)) continue;

			for (size_t i = 0; i < kv.second.size(); i++) {
				uint64_t next = kv.second[i].to;
				uint64_t ek = edge_key(start, next);
				if (visited.count(ek)) continue;
				visited.insert(ek);

				// (K-1)-mer 노드 → 염기 1개 추가로 contig 시작
				string contig = decode_kmer(start, dbg.K - 1);
				contig += decode_base(next);
				uint64_t cur = next;

				// [수정] 직선 구간은 기존대로, 분기점에서는 PE 투표로 경로 선택
				while (true) {
					int id = indeg.count(cur) ? indeg[cur] : 0;
					int od = outdeg.count(cur) ? outdeg[cur] : 0;

					auto it = dbg.adj.find(cur);
					if (it == dbg.adj.end() || it->second.empty()) break;

					if (id == 1 && od == 1) {
						// 분기 없는 직선 구간 — 기존 방식
						uint64_t nxt = it->second[0].to;
						uint64_t nk = edge_key(cur, nxt);
						if (visited.count(nk)) break;
						visited.insert(nk);
						contig += decode_base(nxt);
						cur = nxt;
					}
					else if (od >= 2) {
						// 분기점 — paired-end 투표로 경로 선택
						auto pe_it = pe_links.find(cur);
						if (pe_it == pe_links.end()) break;  // 링크 없으면 끊음

						// 각 후보 경로에 대해 paired-end 지지 수 합산
						uint64_t best_nxt = 0;
						int best_score = 0;
						for (size_t j = 0; j < it->second.size(); j++) {
							uint64_t candidate = it->second[j].to;
							auto score_it = pe_it->second.find(candidate);
							int score = (score_it != pe_it->second.end()) ? score_it->second : 0;
							if (score > best_score) {
								best_score = score;
								best_nxt = candidate;
							}
						}

						// 지지 수 2 미만이면 신뢰도 낮으므로 끊음
						if (best_score < 2 || best_nxt == 0) break;

						uint64_t nk = edge_key(cur, best_nxt);
						if (visited.count(nk)) break;
						visited.insert(nk);
						contig += decode_base(best_nxt);
						cur = best_nxt;
					}
					else {
						break;
					}
				}

				if ((int)contig.size() >= min_len) contigs.push_back(contig);
			}
		}
	}

	// 긴 contig가 앞에 오도록 정렬
	sort(contigs.begin(), contigs.end(),
		[](const string& a, const string& b) { return a.size() > b.size(); });

	return contigs;
}

// contig에 대한 k-mer 인덱스 구축 — 리드 매핑에 사용
struct KmerIdx { int cid, pos; };  // contig id + 위치
struct KmerHit { int cid, pos; };

HashMap<uint64_t, vector<KmerIdx>> build_kmer_index(const vector<string>& contigs, int klen)
{
	HashMap<uint64_t, vector<KmerIdx>> idx;
	idx.reserve(contigs.size() * 200);
	uint64_t mask = (klen < 32) ? (1ULL << (2 * klen)) - 1ULL : ~0ULL;

	for (int ci = 0; ci < (int)contigs.size(); ci++) {
		const string& c = contigs[ci];
		if ((int)c.size() < klen) continue;
		uint64_t val = 0;
		for (int i = 0; i < klen - 1; i++) val = ((val << 2) & mask) | encode_base(c[i]);
		for (int i = klen - 1; i < (int)c.size(); i++) {
			val = ((val << 2) & mask) | encode_base(c[i]);
			KmerIdx ki; ki.cid = ci; ki.pos = i - klen + 1;
			idx[val].push_back(ki);
		}
	}
	return idx;
}

// 리드의 첫 k-mer로 contig 인덱스에서 빠르게 매핑 위치를 탐색
KmerHit map_read(const string& read, const HashMap<uint64_t, vector<KmerIdx>>& idx, int klen)
{
	KmerHit fail; fail.cid = -1; fail.pos = -1;
	if ((int)read.size() < klen) return fail;

	uint64_t mask = (klen < 32) ? (1ULL << (2 * klen)) - 1ULL : ~0ULL;
	uint64_t val = 0;
	for (int i = 0; i < klen; i++) val = ((val << 2) & mask) | encode_base(read[i]);

	auto it = idx.find(val);
	if (it != idx.end() && !it->second.empty()) {
		KmerHit h; h.cid = it->second[0].cid; h.pos = it->second[0].pos;
		return h;
	}
	return fail;
}

// 서로 다른 contig에 매핑된 paired-end 링크 탐지
// insert size로 두 contig 간 gap을 추정하고 지지 수가 충분한 링크만 보고
struct PELink { int from, to, gap_estimate, support; };

vector<PELink> detect_paired_links(
	const vector<PairedRead>& pairs, const vector<string>& contigs,
	const HashMap<uint64_t, vector<KmerIdx>>& idx, int klen, int read_len)
{
	map<pair<int, int>, pair<long long, int>> acc;  // (from,to) → (gap합, 횟수)

	for (size_t i = 0; i < pairs.size(); i++) {
		string read2_forward = reverse_complement(pairs[i].read2_rc);
		KmerHit h1 = map_read(pairs[i].read1, idx, klen);
		KmerHit h2 = map_read(read2_forward, idx, klen);

		// 두 리드가 서로 다른 contig에 매핑된 경우만 유효한 링크
		if (h1.cid < 0 || h2.cid < 0 || h1.cid == h2.cid) continue;

		// insert size - 이미 contig에 포함된 거리로 gap 추정
		int left_remaining = max(0, (int)contigs[h1.cid].size() - (h1.pos + read_len));
		int gap = pairs[i].insert_size - (2 * read_len) - left_remaining - h2.pos;
		gap = max(gap, -(klen - 1));
		gap = min(gap, pairs[i].insert_size);

		auto key = make_pair(h1.cid, h2.cid);
		acc[key].first += gap;
		acc[key].second += 1;
	}

	vector<PELink> links;
	for (auto& kv : acc) {
		int support = kv.second.second;
		if (support < 2) continue;  // 지지 리드 2개 미만 제외
		PELink lk;
		lk.from = kv.first.first; lk.to = kv.first.second;
		lk.gap_estimate = (int)(kv.second.first / support);
		lk.support = support;
		links.push_back(lk);
	}

	// 지지 수 내림차순 정렬
	sort(links.begin(), links.end(),
		[](const PELink& a, const PELink& b) { return a.support > b.support; });
	return links;
}

// scaffold 조립 결과
struct ScaffoldResult {
	string sequence;
	int contig_count_used;
};

// paired-end 링크를 따라 contig들을 이어 scaffold 조립
// gap > 0: 'N' 삽입 / gap < 0: overlap만큼 끝 트리밍
ScaffoldResult build_scaffold(
	const vector<string>& contigs, const vector<PELink>& links,
	int K, int min_support)
{
	ScaffoldResult res = { "", 0 };
	if (contigs.empty()) return res;

	int n = (int)contigs.size();
	vector<bool> used(n, false);
	vector<int> chain, gap_sizes;

	// 각 노드의 최선 입력/출력 링크 선택 (지지 수 최대)
	vector<int> best_out(n, -1), best_out_gap(n, 0), best_out_support(n, -1);
	vector<int> best_in(n, -1);

	for (size_t i = 0; i < links.size(); i++) {
		const PELink& lk = links[i];
		if (lk.support < min_support) continue;
		if (lk.from < 0 || lk.from >= n || lk.to < 0 || lk.to >= n) continue;

		if (lk.support > best_out_support[lk.from]) {
			best_out[lk.from] = lk.to; best_out_gap[lk.from] = lk.gap_estimate;
			best_out_support[lk.from] = lk.support;
		}
		if (best_in[lk.to] == -1) {
			best_in[lk.to] = lk.from;
		}
	}

	// 가장 긴 contig를 앵커로 잡고, 입력 링크를 역추적해 체인 시작점 탐색
	int anchor = 0;
	for (int i = 1; i < n; i++)
		if (contigs[i].size() > contigs[anchor].size()) anchor = i;

	HashSet<int> seen_back;
	int start = anchor;
	while (best_in[start] != -1 && !seen_back.count(best_in[start])) {
		seen_back.insert(start); start = best_in[start];
	}

	// 시작점에서 출력 링크를 따라 체인 구성
	for (int cur = start; cur != -1 && !used[cur]; ) {
		chain.push_back(cur); used[cur] = true;
		int nxt = best_out[cur];
		if (nxt == -1 || used[nxt]) break;
		gap_sizes.push_back(max(best_out_gap[cur], -(K - 1)));
		cur = nxt;
	}
	if (chain.empty()) { chain.push_back(anchor); used[anchor] = true; }

	// 체인에 따라 contig를 이어 붙임 (gap 처리 포함)
	string seq = contigs[chain[0]];
	for (int i = 1; i < (int)chain.size(); i++) {
		int gap = gap_sizes[i - 1];
		if (gap <= 0) {
			// 음수 gap = overlap → 끝 부분 트리밍
			int trim_val = min(-gap, K - 1);
			trim_val = min(trim_val, (int)seq.size());
			seq.resize(seq.size() - trim_val);
		}
		else {
			seq += string(gap, 'N');  // 양수 gap → 'N'으로 채움
		}
		seq += contigs[chain[i]];
	}

	res.sequence = seq; res.contig_count_used = (int)chain.size();
	return res;
}

// ── 정확도 평가 ──────────────────────────────────────────────────────

struct LCSStringMetrics { double accuracy; int matched, original_len, assembly_len; bool reverse_used; };
struct KmerMetrics { double recall, precision, f1; int orig_kmers, asm_kmers, common; };

// LCS(최장 공통 부분수열) 기반 정확도 평가 — 갭/삽입결실 허용
LCSStringMetrics evaluate_lcs_string_accuracy(const string& original, const string& assembly)
{
	auto lcs_one_direction = [&](const string& asm_seq) -> int {
		int n = (int)original.size(), m = (int)asm_seq.size();
		vector<int> prev(m + 1, 0), curr(m + 1, 0);
		for (int i = 1; i <= n; i++) {
			fill(curr.begin(), curr.end(), 0);
			for (int j = 1; j <= m; j++) {
				char lc_a = original[i - 1], lc_b = asm_seq[j - 1];
				if (lc_b == 'N') curr[j] = max(prev[j], curr[j - 1]);
				else if (lc_a == lc_b) curr[j] = prev[j - 1] + 1;
				else                   curr[j] = max(prev[j], curr[j - 1]);
			}
			prev.swap(curr);
		}
		return prev[m];
	};

	int fwd = lcs_one_direction(assembly);
	int rev = lcs_one_direction(reverse_complement(assembly));
	int lc_matched = max(fwd, rev);

	LCSStringMetrics result;
	result.accuracy = ((int)original.size() > 0)
		? (double)lc_matched / original.size() * 100.0 : 0.0;
	result.matched = lc_matched; result.original_len = (int)original.size();
	result.assembly_len = (int)assembly.size(); result.reverse_used = (rev > fwd);
	return result;
}

// k-mer 집합 비교로 recall / precision / F1 계산
KmerMetrics evaluate_kmer(const string& original, const string& assembly, int k)
{
	auto make_kmer_set = [&](const string& seq) -> HashSet<uint64_t> {
		HashSet<uint64_t> s;
		int n = (int)seq.size();
		if (n < k) return s;
		uint64_t mask = (k < 32) ? (1ULL << (2 * k)) - 1ULL : ~0ULL;
		uint64_t val = 0; int valid_count = 0;
		for (int i = 0; i < n; i++) {
			char c = seq[i];
			if (c == 'N') { val = 0; valid_count = 0; continue; }
			val = ((val << 2) & mask) | encode_base(c);
			if (++valid_count >= k) s.insert(val);
		}
		return s;
	};

	HashSet<uint64_t> orig = make_kmer_set(original);
	HashSet<uint64_t> asm_set = make_kmer_set(assembly);
	int lc_common = 0;
	for (const auto& x : asm_set) if (orig.count(x)) lc_common++;

	KmerMetrics km;
	km.orig_kmers = (int)orig.size(); km.asm_kmers = (int)asm_set.size(); km.common = lc_common;
	km.recall = orig.empty() ? 0.0 : (double)lc_common / orig.size() * 100.0;
	km.precision = asm_set.empty() ? 0.0 : (double)lc_common / asm_set.size() * 100.0;
	km.f1 = (km.recall + km.precision > 0.0)
		? 2.0 * km.recall * km.precision / (km.recall + km.precision) : 0.0;
	return km;
}

// contig 길이 기준 N50 계산 (전체 길이의 50%를 커버하는 최소 contig 길이)
int calculate_n50(const vector<string>& contigs)
{
	if (contigs.empty()) return 0;
	vector<int> lens; int total = 0;
	for (size_t i = 0; i < contigs.size(); i++) { lens.push_back((int)contigs[i].size()); total += lens.back(); }
	sort(lens.rbegin(), lens.rend());
	int sum = 0;
	for (size_t i = 0; i < lens.size(); i++) { sum += lens[i]; if (sum * 2 >= total) return lens[i]; }
	return 0;
}

long long elapsed_ms(
	chrono::high_resolution_clock::time_point start,
	chrono::high_resolution_clock::time_point end)
{
	return chrono::duration_cast<chrono::milliseconds>(end - start).count();
}

int main(int argc, char* argv[])
{
	init_tables();
	srand((unsigned)time(0));
	cout << fixed << setprecision(4);

	// ── 주요 파라미터 ──────────────────────────────────
	const int K = 21;
	const int PAIR_COUNT = 50000;
	const int READ_LEN = 100;
	const int INSERT_MIN = 200;
	const int INSERT_MAX = 1000;
	const int MAX_MISMATCH = 0;
	const int MIN_WEIGHT = 2;
	const int MIN_CONTIG = 50;
	const int MIN_PE_SUPPORT = 2;
	const size_t GENOME_LEN = 100000;

	string genome_path = (argc >= 2) ? argv[1] : "../chr22.fa";
	auto load_start = chrono::high_resolution_clock::now();
	string genome = load_fasta(genome_path, GENOME_LEN);
	auto load_end = chrono::high_resolution_clock::now();
	long long runtime_load = elapsed_ms(load_start, load_end);
	if (genome.empty()) { cerr << "genome 로드 실패: " << genome_path << "\n"; return 1; }

	// ── 실험 환경 출력 ─────────────────────────────────
	cout << "========================================\n";
	cout << "실험 환경\n";
	cout << "========================================\n";
	cout << "genome 길이    : " << genome.size() << " bp\n";
	cout << "reads 수       : " << PAIR_COUNT * 2 << "개 (paired-end " << PAIR_COUNT << "쌍)\n";
	cout << "read 길이      : " << READ_LEN << " bp\n";
	cout << "mismatch 허용  : " << MAX_MISMATCH << "개\n";
	cout << "insert size 범위: " << INSERT_MIN << " ~ " << INSERT_MAX << " bp\n";
	cout << "\n";

	auto assembly_start = chrono::high_resolution_clock::now();

	// 1. paired-end 리드 시뮬레이션
	vector<PairedRead> paired_reads = generate_paired_reads(
		genome, PAIR_COUNT, READ_LEN, INSERT_MIN, INSERT_MAX, MAX_MISMATCH);

	// 2. De Bruijn 그래프 구축 → 저빈도 간선 가지치기
	DeBruijnGraph dbg(K);
	for (size_t i = 0; i < paired_reads.size(); i++) {
		dbg.add_read(paired_reads[i].read1);
		dbg.add_read(reverse_complement(paired_reads[i].read2_rc));
	}
	dbg.prune(MIN_WEIGHT);

	// 3. [추가] paired-end 링크 테이블 구축 → contig 추출 (PE guided)
	PELinkTable pe_link_table = build_pe_link_table(paired_reads, K);
	vector<string> contigs = generate_contigs_unitig(dbg, MIN_CONTIG, pe_link_table);

	// 4. contig 인덱스 구축 → paired-end 링크 탐지 → scaffold 조립
	HashMap<uint64_t, vector<KmerIdx>> kmer_idx = build_kmer_index(contigs, K);
	vector<PELink> links = detect_paired_links(paired_reads, contigs, kmer_idx, K, READ_LEN);
	ScaffoldResult scaffold = build_scaffold(contigs, links, K, MIN_PE_SUPPORT);

	auto assembly_end = chrono::high_resolution_clock::now();
	long long runtime_assembly = elapsed_ms(assembly_start, assembly_end);

	// 5. 품질 평가 (LCS / k-mer 지표)
	auto eval_start = chrono::high_resolution_clock::now();
	LCSStringMetrics lsm = evaluate_lcs_string_accuracy(genome, scaffold.sequence);
	KmerMetrics      km = evaluate_kmer(genome, scaffold.sequence, K);
	auto eval_end = chrono::high_resolution_clock::now();
	long long runtime_eval = elapsed_ms(eval_start, eval_end);

	size_t total_contig = 0, longest = 0;
	for (size_t i = 0; i < contigs.size(); i++) {
		total_contig += contigs[i].size();
		if (contigs[i].size() > longest) longest = contigs[i].size();
	}
	int n50 = calculate_n50(contigs);
	double memory_mb = get_memory_usage_mb();

	// ── 결과 출력 ──────────────────────────────────────
	cout << "\n========================================\n";
	cout << "  De Bruijn Graph Assembly Result\n";
	cout << "========================================\n";

	cout << "\n[Resource]\n";
	cout << "  Memory Usage         : " << memory_mb << " MB\n";

	cout << "\n[Contig]\n";
	cout << "  Contig Count         : " << contigs.size() << "\n";
	cout << "  Total Contig Length  : " << total_contig << " bp\n";
	cout << "  Longest Contig       : " << longest << " bp\n";
	cout << "  N50                  : " << n50 << " bp\n";

	cout << "\n[Accuracy]\n";
	cout << "  정확도               : " << lsm.accuracy << " %\n";
	cout << "  Correct Bases        : " << lsm.matched << " / " << lsm.original_len << "\n";
	cout << "  Missing/Wrong Bases  : " << (lsm.original_len - lsm.matched) << "\n";
	cout << "  k-mer Recall         : " << km.recall << " %\n";
	cout << "  k-mer Precision      : " << km.precision << " %\n";
	cout << "  k-mer F1             : " << km.f1 << " %\n";

	{
		const int BAR_LEN = 25;
		auto print_bar = [&](const string& label, double value, double max_val) {
			int filled = (max_val > 0.0) ? (int)(value / max_val * BAR_LEN) : 0;
			if (filled > BAR_LEN) filled = BAR_LEN;
			cout << "  " << left << setw(20) << label << " |";
			for (int i = 0; i < BAR_LEN; i++) cout << (i < filled ? "#" : " ");
			cout << "| " << fixed << setprecision(1) << value << "%\n";
		};

		cout << "\n[정확도 시각화]\n";
		cout << "  " << string(48, '-') << "\n";
		cout << "  " << left << setw(20) << "" << "  0%    25%   50%   75%  100%\n";
		print_bar("정확도", lsm.accuracy, 100.0);
		print_bar("k-mer Recall", km.recall, 100.0);
		print_bar("k-mer Precision", km.precision, 100.0);
		print_bar("k-mer F1", km.f1, 100.0);
		cout << "  " << string(48, '-') << "\n";
	}

	cout << "\n[실행시간]\n";
	cout << "  Load Runtime         : " << runtime_load << " ms\n";
	cout << "  Assembly Runtime     : " << runtime_assembly << " ms\n";
	cout << "  Eval Runtime         : " << runtime_eval << " ms\n";
	cout << "  Total Runtime        : " << runtime_load + runtime_assembly + runtime_eval << " ms\n";

	cout << "========================================\n\n";

	return 0;
}