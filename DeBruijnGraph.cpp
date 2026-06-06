#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <chrono>
#include <random>
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

static bool DNA_VALID[256];
static char DNA_UPPER[256];
static uint8_t BASE_ENC[256];
static const char BASE_DEC[4] = { 'A', 'C', 'G', 'T' };

static void init_tables()
{
	for (int i = 0; i < 256; i++) {
		DNA_VALID[i] = false;
		DNA_UPPER[i] = (char)i;
		BASE_ENC[i] = 0;
	}

	DNA_VALID['A'] = DNA_VALID['C'] = DNA_VALID['G'] = DNA_VALID['T'] = true;
	DNA_VALID['a'] = DNA_VALID['c'] = DNA_VALID['g'] = DNA_VALID['t'] = true;

	DNA_UPPER['a'] = 'A';
	DNA_UPPER['c'] = 'C';
	DNA_UPPER['g'] = 'G';
	DNA_UPPER['t'] = 'T';

	BASE_ENC['A'] = BASE_ENC['a'] = 0;
	BASE_ENC['C'] = BASE_ENC['c'] = 1;
	BASE_ENC['G'] = BASE_ENC['g'] = 2;
	BASE_ENC['T'] = BASE_ENC['t'] = 3;
}

double get_memory_usage_mb()
{
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS pmc;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
		return (double)pmc.WorkingSetSize / (1024.0 * 1024.0);
	}
	return 0.0;
#else
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage) == 0) {
		return (double)usage.ru_maxrss / 1024.0;
	}
	return 0.0;
#endif
}

inline uint64_t encode_base(char c) { return BASE_ENC[(unsigned char)c]; }
inline char     decode_base(uint64_t x) { return BASE_DEC[x & 3]; }

string decode_kmer(uint64_t value, int length)
{
	string s(length, 'A');
	for (int i = length - 1; i >= 0; i--) { s[i] = decode_base(value); value >>= 2; }
	return s;
}

string load_fasta(const string& filename, size_t max_len = 0)
{
	ifstream f(filename);

	if (!f.is_open()) {
		cerr << "FASTA open fail: " << filename << "\n";
		return "";
	}

	string genome;
	string line;

	if (max_len > 0) {
		genome.reserve(max_len);
	}

	while (getline(f, line)) {
		if (line.empty() || line[0] == '>') {
			continue;
		}

		for (unsigned char c : line) {
			if (DNA_VALID[c]) {
				genome += DNA_UPPER[c];

				if (max_len > 0 && genome.size() >= max_len) {
					genome.shrink_to_fit();
					return genome;
				}
			}
		}
	}

	genome.shrink_to_fit();
	return genome;
}

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

void add_random_mismatch(string& read, int max_mismatch)
{
	if (max_mismatch <= 0 || read.empty()) {
		return;
	}

	int num_mutations = rand() % (max_mismatch + 1);
	string bases = "ACGT";

	for (int m = 0; m < num_mutations; m++) {
		int pos = rand() % (int)read.size();
		char original = read[pos];
		char mutated;

		do {
			mutated = bases[rand() % 4];
		} while (mutated == original);

		read[pos] = mutated;
	}
}

struct PairedRead {
	string read1;
	string read2_rc;
	int frag_start;
	int insert_size;
};

vector<PairedRead> generate_paired_reads(
	const string& genome,
	int pair_count,
	int read_len,
	int insert_min,
	int insert_max,
	int max_mismatch
)
{
	vector<PairedRead> pairs;
	pairs.reserve(pair_count);

	int n = (int)genome.size();

	for (int i = 0; i < pair_count; i++) {
		int frag_len = insert_min + rand() % (insert_max - insert_min);

		if (frag_len + read_len >= n) {
			continue;
		}

		int frag_start = rand() % (n - frag_len - read_len);

		string read1 = genome.substr(frag_start, read_len);
		add_random_mismatch(read1, max_mismatch);

		string read2 = genome.substr(frag_start + frag_len, read_len);
		string read2_rc = reverse_complement(read2);
		add_random_mismatch(read2_rc, max_mismatch);

		PairedRead pr;
		pr.read1 = read1;
		pr.read2_rc = read2_rc;
		pr.frag_start = frag_start;
		pr.insert_size = frag_len + read_len;

		pairs.push_back(pr);
	}

	return pairs;
}

bool encode_kmer_from_string(const string& s, int start, int k, uint64_t& value)
{
	if (start < 0 || start + k >(int)s.size()) return false;

	value = 0;

	for (int i = 0; i < k; i++) {
		unsigned char c = (unsigned char)s[start + i];

		if (!DNA_VALID[c]) return false;

		value = (value << 2) | encode_base((char)c);
	}

	return true;
}

HashMap<uint64_t, int> count_read_kmers(const vector<PairedRead>& pairs, int k)
{
	HashMap<uint64_t, int> counts;
	counts.reserve(pairs.size() * 100);

	for (size_t i = 0; i < pairs.size(); i++) {
		string read2_forward = reverse_complement(pairs[i].read2_rc);
		const string* reads[2] = { &pairs[i].read1, &read2_forward };

		for (int r = 0; r < 2; r++) {
			const string& read = *reads[r];

			if ((int)read.size() < k) continue;

			for (int pos = 0; pos + k <= (int)read.size(); pos++) {
				uint64_t val = 0;

				if (encode_kmer_from_string(read, pos, k, val)) counts[val]++;
			}
		}
	}

	return counts;
}

int local_bad_kmer_score(
	const string& read, int changed_pos, int k,
	const HashMap<uint64_t, int>& kmer_counts, int solid_min)
{
	int n = (int)read.size();

	if (n < k) return 0;

	int start_min = max(0, changed_pos - k + 1);
	int start_max = min(changed_pos, n - k);
	int bad = 0;

	for (int start = start_min; start <= start_max; start++) {
		uint64_t val = 0;

		if (!encode_kmer_from_string(read, start, k, val)) { bad++; continue; }

		auto it = kmer_counts.find(val);

		if (it == kmer_counts.end() || it->second < solid_min) bad++;
	}

	return bad;
}

string correct_read_by_solid_kmers(
	string read, int k,
	const HashMap<uint64_t, int>& kmer_counts,
	int solid_min, int max_correction)
{
	if (max_correction <= 0 || (int)read.size() < k) return read;

	static const char BASES[4] = { 'A', 'C', 'G', 'T' };

	for (int edit = 0; edit < max_correction; edit++) {
		int best_pos = -1;
		char best_base = 0;
		int best_gain = 0;

		for (int pos = 0; pos < (int)read.size(); pos++) {
			char original = read[pos];

			if (!DNA_VALID[(unsigned char)original]) continue;

			int before = local_bad_kmer_score(read, pos, k, kmer_counts, solid_min);

			for (int b = 0; b < 4; b++) {
				char candidate = BASES[b];

				if (candidate == original) continue;

				read[pos] = candidate;

				int after = local_bad_kmer_score(read, pos, k, kmer_counts, solid_min);
				int gain = before - after;

				if (gain > best_gain) { best_gain = gain; best_pos = pos; best_base = candidate; }
			}

			read[pos] = original;
		}

		if (best_pos < 0 || best_gain <= 0) break;

		read[best_pos] = best_base;
	}

	return read;
}

void correct_paired_reads(
	vector<PairedRead>& pairs, int k,
	const HashMap<uint64_t, int>& kmer_counts,
	int solid_min, int max_correction)
{
	for (size_t i = 0; i < pairs.size(); i++) {
		string old_read2_forward = reverse_complement(pairs[i].read2_rc);

		pairs[i].read1 = correct_read_by_solid_kmers(pairs[i].read1, k, kmer_counts, solid_min, max_correction);
		pairs[i].read2_rc = reverse_complement(
			correct_read_by_solid_kmers(old_read2_forward, k, kmer_counts, solid_min, max_correction));
	}
}

struct Edge {
	uint64_t to;
	uint32_t weight;
};

class DeBruijnGraph
{
public:
	HashMap<uint64_t, vector<Edge>> adj;
	int K;
	uint64_t node_mask;

	explicit DeBruijnGraph(int k) : K(k)
	{
		assert(k >= 2 && k <= 32);

		if (K - 1 < 32) node_mask = (1ULL << (2 * (K - 1))) - 1ULL;
		else             node_mask = ~0ULL;

		adj.reserve(1 << 20);
	}

	void add_read(const string& read)
	{
		int n = (int)read.size();

		if (n < K) return;

		uint64_t prefix = 0;

		for (int i = 0; i < K - 1; i++)
			prefix = ((prefix << 2) & node_mask) | encode_base(read[i]);

		if (adj.find(prefix) == adj.end()) adj[prefix] = vector<Edge>();

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

vector<string> generate_contigs_unitig(DeBruijnGraph& dbg, int min_len)
{
	HashMap<uint64_t, int> indeg, outdeg;

	for (auto& kv : dbg.adj) {
		uint64_t from = kv.first;
		outdeg[from] = (int)kv.second.size();
		if (indeg.find(from) == indeg.end()) indeg[from] = 0;
		for (size_t i = 0; i < kv.second.size(); i++) indeg[kv.second[i].to]++;
	}

	auto edge_key = [](uint64_t from, uint64_t to) -> uint64_t {
		return from ^ (to * 2654435761ULL);
	};

	auto is_start = [&](uint64_t v) -> bool {
		int id = indeg.count(v) ? indeg[v] : 0;
		int od = outdeg.count(v) ? outdeg[v] : 0;
		if (od == 0) return false;
		return !(id == 1 && od == 1);
	};

	HashSet<uint64_t> visited;
	vector<string> contigs;

	for (int pass = 0; pass < 2; pass++) {
		for (auto& kv : dbg.adj) {
			uint64_t start = kv.first;
			if (pass == 0 && !is_start(start)) continue;

			for (size_t i = 0; i < kv.second.size(); i++) {
				uint64_t next = kv.second[i].to;
				uint64_t ek = edge_key(start, next);

				if (visited.count(ek)) continue;
				visited.insert(ek);

				string contig = decode_kmer(start, dbg.K - 1);
				contig += decode_base(next);
				uint64_t cur = next;

				while (true) {
					int id = indeg.count(cur) ? indeg[cur] : 0;
					int od = outdeg.count(cur) ? outdeg[cur] : 0;
					if (!(id == 1 && od == 1)) break;

					auto it = dbg.adj.find(cur);
					if (it == dbg.adj.end() || it->second.empty()) break;

					uint64_t nxt = it->second[0].to;
					uint64_t nk = edge_key(cur, nxt);
					if (visited.count(nk)) break;

					visited.insert(nk);
					contig += decode_base(nxt);
					cur = nxt;
				}

				if ((int)contig.size() >= min_len) contigs.push_back(contig);
			}
		}
	}

	sort(contigs.begin(), contigs.end(),
		[](const string& a, const string& b) { return a.size() > b.size(); });

	return contigs;
}

struct KmerIdx { int cid, pos; };
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

struct PELink {
	int from, to;
	int gap_estimate;
	int support;
};

vector<PELink> detect_paired_links(
	const vector<PairedRead>& pairs,
	const vector<string>& contigs,
	const HashMap<uint64_t, vector<KmerIdx>>& idx,
	int klen, int read_len)
{
	map<pair<int, int>, pair<long long, int>> acc;

	for (size_t i = 0; i < pairs.size(); i++) {
		string read2_forward = reverse_complement(pairs[i].read2_rc);

		KmerHit h1 = map_read(pairs[i].read1, idx, klen);
		KmerHit h2 = map_read(read2_forward, idx, klen);

		if (h1.cid < 0 || h2.cid < 0 || h1.cid == h2.cid) continue;

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
		if (support < 2) continue;

		PELink lk;
		lk.from = kv.first.first;
		lk.to = kv.first.second;
		lk.gap_estimate = (int)(kv.second.first / support);
		lk.support = support;
		links.push_back(lk);
	}

	sort(links.begin(), links.end(),
		[](const PELink& a, const PELink& b) { return a.support > b.support; });

	return links;
}

struct ScaffoldResult {
	string sequence;
	int contig_count_used, unused_contigs;
	int gap_count, total_gap_length, total_overlap_trimmed;
};

ScaffoldResult build_scaffold(
	const vector<string>& contigs, const vector<PELink>& links,
	int K, int min_support)
{
	ScaffoldResult res;
	res.sequence = "";
	res.contig_count_used = 0;
	res.unused_contigs = 0;
	res.gap_count = 0;
	res.total_gap_length = 0;
	res.total_overlap_trimmed = 0;

	if (contigs.empty()) return res;

	int n = (int)contigs.size();

	vector<bool> used(n, false);
	vector<int>  chain, gap_sizes;
	vector<int>  best_out(n, -1), best_out_gap(n, 0), best_out_support(n, -1);
	vector<int>  best_in(n, -1), best_in_support(n, -1);

	for (size_t i = 0; i < links.size(); i++) {
		const PELink& lk = links[i];
		if (lk.support < min_support) continue;
		if (lk.from < 0 || lk.from >= n || lk.to < 0 || lk.to >= n) continue;

		if (lk.support > best_out_support[lk.from]) {
			best_out[lk.from] = lk.to; best_out_gap[lk.from] = lk.gap_estimate;
			best_out_support[lk.from] = lk.support;
		}
		if (lk.support > best_in_support[lk.to]) {
			best_in[lk.to] = lk.from; best_in_support[lk.to] = lk.support;
		}
	}

	int anchor = 0;
	for (int i = 1; i < n; i++)
		if (contigs[i].size() > contigs[anchor].size()) anchor = i;

	HashSet<int> seen_back;
	int start = anchor;
	while (best_in[start] != -1 && !seen_back.count(best_in[start])) {
		seen_back.insert(start);
		start = best_in[start];
	}

	for (int cur = start; cur != -1 && !used[cur]; ) {
		chain.push_back(cur);
		used[cur] = true;
		int nxt = best_out[cur];
		if (nxt == -1 || used[nxt]) break;
		gap_sizes.push_back(max(best_out_gap[cur], -(K - 1)));
		cur = nxt;
	}

	if (chain.empty()) { chain.push_back(anchor); used[anchor] = true; }

	string seq = contigs[chain[0]];

	for (int i = 1; i < (int)chain.size(); i++) {
		int gap = gap_sizes[i - 1];

		if (gap <= 0) {
			int trim_val = -gap;
			if (trim_val > K - 1) trim_val = K - 1;
			if (trim_val > (int)seq.size()) trim_val = (int)seq.size();
			seq.resize(seq.size() - trim_val);
			res.total_overlap_trimmed += trim_val;
		}
		else {
			seq += string(gap, 'N');
			res.gap_count++;
			res.total_gap_length += gap;
		}

		seq += contigs[chain[i]];
	}

	int unused = 0;
	for (int i = 0; i < n; i++) if (!used[i]) unused++;

	res.sequence = seq;
	res.contig_count_used = (int)chain.size();
	res.unused_contigs = unused;

	return res;
}

// ── 정확도 평가 함수들 ──────────────────────────────────────────────

struct BestStringMetrics { double accuracy; int matched, mismatches, n_bases, compared, original_len, assembly_len, best_offset; bool reverse_used; };
struct LCSStringMetrics { double accuracy; int matched, original_len, assembly_len; bool reverse_used; };
struct KmerMetrics { double recall, precision, f1; int orig_kmers, asm_kmers, common; };

BestStringMetrics evaluate_best_string_alignment(const string& original, const string& assembly)
{
	BestStringMetrics best;
	best.accuracy = 0.0;
	best.matched = 0;
	best.mismatches = 0;
	best.n_bases = 0;
	best.compared = 0;
	best.original_len = (int)original.size();
	best.assembly_len = (int)assembly.size();
	best.best_offset = 0;
	best.reverse_used = false;

	// 람다 내부 변수명을 외부 변수와 구분되도록 lc_ 접두사 사용
	auto test_one_direction = [&](const string& asm_seq, bool reverse_flag) {
		int olen = (int)original.size(), alen = (int)asm_seq.size();

		for (int offset = -alen + 1; offset < olen; offset++) {
			int lc_matched = 0, lc_mismatches = 0, lc_n_bases = 0, lc_compared = 0;

			for (int ai = 0; ai < alen; ai++) {
				int oi = offset + ai;
				if (oi < 0 || oi >= olen) continue;
				char a = asm_seq[ai];
				if (a == 'N') { lc_n_bases++; continue; }
				lc_compared++;
				if (original[oi] == a) lc_matched++; else lc_mismatches++;
			}

			if (lc_compared == 0) continue;

			double acc = (double)lc_matched / (double)olen * 100.0;
			if (acc > best.accuracy) {
				best.accuracy = acc;
				best.matched = lc_matched;
				best.mismatches = lc_mismatches;
				best.n_bases = lc_n_bases;
				best.compared = lc_compared;
				best.original_len = (int)original.size();
				best.assembly_len = (int)asm_seq.size();
				best.best_offset = offset;
				best.reverse_used = reverse_flag;
			}
		}
	};

	test_one_direction(assembly, false);
	test_one_direction(reverse_complement(assembly), true);

	return best;
}

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
				else               curr[j] = max(prev[j], curr[j - 1]);
			}
			prev.swap(curr);
		}

		return prev[m];
	};

	int fwd = lcs_one_direction(assembly);
	int rev = lcs_one_direction(reverse_complement(assembly));
	int lc_matched = max(fwd, rev);

	LCSStringMetrics result;
	result.accuracy = ((int)original.size() > 0) ? (double)lc_matched / original.size() * 100.0 : 0.0;
	result.matched = lc_matched;
	result.original_len = (int)original.size();
	result.assembly_len = (int)assembly.size();
	result.reverse_used = (rev > fwd);
	return result;
}

KmerMetrics evaluate_kmer(const string& original, const string& assembly, int k)
{
	auto make_kmer_set = [&](const string& seq) -> HashSet<uint64_t> {
		HashSet<uint64_t> s;
		int n = (int)seq.size();
		if (n < k) return s;
		uint64_t mask = (k < 32) ? (1ULL << (2 * k)) - 1ULL : ~0ULL;
		uint64_t val = 0;
		int valid_count = 0;
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
	km.orig_kmers = (int)orig.size();
	km.asm_kmers = (int)asm_set.size();
	km.common = lc_common;
	km.recall = orig.empty() ? 0.0 : (double)lc_common / orig.size() * 100.0;
	km.precision = asm_set.empty() ? 0.0 : (double)lc_common / asm_set.size() * 100.0;
	km.f1 = (km.recall + km.precision > 0.0)
		? 2.0 * km.recall * km.precision / (km.recall + km.precision) : 0.0;

	return km;
}

int calculate_n50(const vector<string>& contigs)
{
	if (contigs.empty()) return 0;

	vector<int> lens;
	int total = 0;
	for (size_t i = 0; i < contigs.size(); i++) { lens.push_back((int)contigs[i].size()); total += lens.back(); }
	sort(lens.rbegin(), lens.rend());

	int sum = 0;
	for (size_t i = 0; i < lens.size(); i++) {
		sum += lens[i];
		if (sum * 2 >= total) return lens[i];
	}

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

	const int K = 21;
	const int PAIR_COUNT = 30000;
	const int READ_LEN = 100;
	const int INSERT_MIN = 200;
	const int INSERT_MAX = 400;
	const int MAX_MISMATCH = 1;
	const int MIN_WEIGHT = 10;
	const int MIN_CONTIG = 200;
	const int MIN_PE_SUPPORT = 2;
	const size_t GENOME_LEN = 10000;

	string genome_path = (argc >= 2) ? argv[1] : "../chr22.fa";
	string genome = load_fasta(genome_path, GENOME_LEN);

	if (genome.empty()) { cerr << "genome 로드 실패: " << genome_path << "\n"; return 1; }

	auto assembly_start = chrono::high_resolution_clock::now();

	// 1. paired-end read 시뮬레이션
	vector<PairedRead> paired_reads = generate_paired_reads(
		genome, PAIR_COUNT, READ_LEN, INSERT_MIN, INSERT_MAX, MAX_MISMATCH);

	// 2. De Bruijn 그래프 구축 → 가지치기 → contig 추출
	DeBruijnGraph dbg(K);
	for (size_t i = 0; i < paired_reads.size(); i++) {
		dbg.add_read(paired_reads[i].read1);
		dbg.add_read(reverse_complement(paired_reads[i].read2_rc));
	}
	dbg.prune(MIN_WEIGHT);
	vector<string> contigs = generate_contigs_unitig(dbg, MIN_CONTIG);

	// 3. paired-end 링크 탐지 → scaffold 조립
	HashMap<uint64_t, vector<KmerIdx>> kmer_idx = build_kmer_index(contigs, K);
	vector<PELink> links = detect_paired_links(paired_reads, contigs, kmer_idx, K, READ_LEN);
	ScaffoldResult scaffold = build_scaffold(contigs, links, K, MIN_PE_SUPPORT);

	auto assembly_end = chrono::high_resolution_clock::now();
	long long runtime_assembly = elapsed_ms(assembly_start, assembly_end);

	// 4. 품질 평가
	BestStringMetrics bsm = evaluate_best_string_alignment(genome, scaffold.sequence);
	LCSStringMetrics  lsm = evaluate_lcs_string_accuracy(genome, scaffold.sequence);
	KmerMetrics       km = evaluate_kmer(genome, scaffold.sequence, K);

	size_t total_contig = 0, longest = 0;
	for (size_t i = 0; i < contigs.size(); i++) {
		total_contig += contigs[i].size();
		if (contigs[i].size() > longest) longest = contigs[i].size();
	}

	int n50 = calculate_n50(contigs);
	double memory_mb = get_memory_usage_mb();

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
	cout << "  LCS String Accuracy  : " << lsm.accuracy << " %\n";
	cout << "  Correct Bases        : " << lsm.matched << " / " << lsm.original_len << "\n";
	cout << "  Missing/Wrong Bases  : " << (lsm.original_len - lsm.matched) << "\n";
	cout << "  Best String Accuracy : " << bsm.accuracy << " %\n";
	cout << "  Best Offset          : " << bsm.best_offset << "\n";
	cout << "  k-mer Recall         : " << km.recall << " %\n";
	cout << "  k-mer Precision      : " << km.precision << " %\n";
	cout << "  k-mer F1             : " << km.f1 << " %\n";

	// ── 정확도 시각화 ──────────────────────────────────────────────
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
	cout << "  Assembly Runtime     : " << runtime_assembly << " ms\n";

	cout << "========================================\n\n";

	return 0;
}