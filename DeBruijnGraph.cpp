#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cctype>
#include <chrono>
#include <random>
#include <unordered_map>
#include <unordered_set>

template<class K, class V> using HashMap = std::unordered_map<K, V>;
template<class K>          using HashSet = std::unordered_set<K>;

using namespace std;

// ──────────────────────────────────────────
// 공통 lookup 테이블
// ──────────────────────────────────────────
static bool     DNA_VALID[256];
static char     DNA_UPPER[256];
static uint8_t  BASE_ENC[256];
static const char BASE_DEC[4] = { 'A','C','G','T' };

static void init_tables()
{
	for (int i = 0; i < 256; i++) { DNA_VALID[i] = false; DNA_UPPER[i] = (char)i; }
	DNA_VALID['A'] = DNA_VALID['C'] = DNA_VALID['G'] = DNA_VALID['T'] = true;
	DNA_VALID['a'] = DNA_VALID['c'] = DNA_VALID['g'] = DNA_VALID['t'] = true;
	DNA_UPPER['a'] = 'A'; DNA_UPPER['c'] = 'C';
	DNA_UPPER['g'] = 'G'; DNA_UPPER['t'] = 'T';

	for (int i = 0; i < 256; i++) BASE_ENC[i] = 0;
	BASE_ENC['A'] = 0; BASE_ENC['C'] = 1; BASE_ENC['G'] = 2; BASE_ENC['T'] = 3;
}

inline uint64_t encode_base(char c) { return BASE_ENC[(unsigned char)c]; }
inline char     decode_base(uint64_t x) { return BASE_DEC[x & 3]; }

string decode_kmer(uint64_t value, int length)
{
	string result(length, 'A');
	for (int i = length - 1; i >= 0; i--)
	{
		result[i] = decode_base(value);
		value >>= 2;
	}
	return result;
}

// ──────────────────────────────────────────
// 1. FASTA 로드
// ──────────────────────────────────────────
string load_fasta(const string& filename)
{
	ifstream file(filename);
	if (!file.is_open()) { cerr << "FASTA 파일 열기 실패: " << filename << "\n"; return ""; }

	string line, genome;
	while (getline(file, line))
	{
		if (line.empty() || line[0] == '>') continue;
		for (unsigned char c : line)
			if (DNA_VALID[c]) genome += DNA_UPPER[c];
	}
	genome.shrink_to_fit();
	return genome;
}

// ──────────────────────────────────────────
// 2. De Bruijn Graph
//    uint64_t 인코딩으로 k <= 32 지원
// ──────────────────────────────────────────
struct Edge { uint64_t suffix; uint32_t weight; };  // suffix: uint64_t로 변경

class DeBruijnGraph
{
private:
	HashMap<uint64_t, vector<Edge>> graph;  // key: uint64_t로 변경
	int k;
	uint64_t mask;  // uint64_t로 변경

public:
	DeBruijnGraph(int kmer_size) : k(kmer_size)
	{
		if (k > 32) { cerr << "이 버전은 k <= 32\n"; exit(1); }  // 상한 32로 확대
		if (k < 2) { cerr << "k >= 2 이어야 합니다\n"; exit(1); }
		// 1ULL: uint32_t 오버플로 방지
		mask = (1ULL << (2 * (k - 1))) - 1ULL;
		graph.reserve(1 << 22);  // k 커질수록 고유 노드 증가 → 예약 크기 확대
	}

	void add_read(const string& read)
	{
		const int n = (int)read.length();
		if (n < k) return;

		uint64_t prefix = 0;  // uint64_t로 변경
		for (int i = 0; i < k - 1; i++)
		{
			prefix <<= 2;
			prefix |= encode_base(read[i]);
		}
		for (int i = k - 1; i < n; i++)
		{
			uint64_t suffix = ((prefix << 2) & mask)  // uint64_t로 변경
				| encode_base(read[i]);

			auto& edges = graph[prefix];
			bool found = false;
			for (auto& e : edges)
				if (e.suffix == suffix) { e.weight++; found = true; break; }
			if (!found) edges.push_back({ suffix, 1 });

			if (graph.find(suffix) == graph.end())
				graph.emplace(suffix, vector<Edge>{});

			prefix = suffix;
		}
	}

	void prune(int min_weight = 2)
	{
		for (auto& kv : graph)
		{
			auto& edges = kv.second;
			edges.erase(
				remove_if(edges.begin(), edges.end(),
					[min_weight](const Edge& e) { return (int)e.weight < min_weight; }),
				edges.end());
		}
	}

	void build_degree_cache(
		HashMap<uint64_t, int>& indeg,   // uint64_t로 변경
		HashMap<uint64_t, int>& outdeg) const
	{
		indeg.reserve(graph.size());
		outdeg.reserve(graph.size());
		for (const auto& kv : graph)
		{
			outdeg[kv.first] = (int)kv.second.size();
			for (const auto& e : kv.second) indeg[e.suffix]++;
		}
	}

	uint64_t get_best_next(uint64_t node) const  // uint64_t로 변경
	{
		const auto& edges = graph.at(node);
		uint64_t best = 0, best_w = 0;  // uint64_t로 변경
		for (const auto& e : edges)
			if (e.weight > best_w) { best_w = e.weight; best = e.suffix; }
		return best;
	}

	vector<string> generate_contigs(int min_len = 0)
	{
		HashMap<uint64_t, int> indeg, outdeg;  // uint64_t로 변경
		build_degree_cache(indeg, outdeg);

		vector<string> contigs;
		contigs.reserve(4096);

		for (const auto& kv : graph)
		{
			uint64_t start = kv.first;  // uint64_t로 변경
			if (indeg[start] == 1 && outdeg[start] == 1) continue;

			for (const auto& first_edge : kv.second)
			{
				uint64_t current = first_edge.suffix;  // uint64_t로 변경
				string contig = decode_kmer(start, k - 1);
				contig.reserve(512);
				contig += decode_base(current);

				while (indeg[current] == 1 && outdeg[current] == 1)
				{
					current = get_best_next(current);
					contig += decode_base(current);
				}
				if ((int)contig.length() >= min_len)
					contigs.push_back(contig);
			}
		}

		if (contigs.empty() && !graph.empty())
		{
			uint64_t start = graph.begin()->first;  // uint64_t로 변경
			uint64_t current = start;
			string contig = decode_kmer(start, k - 1);
			do { current = get_best_next(current); contig += decode_base(current); } while (current != start);
			contigs.push_back(contig);
		}

		return contigs;
	}
};

// ──────────────────────────────────────────
// 3. FASTQ 스트리밍
// ──────────────────────────────────────────
long long stream_fastq(const string& filename, DeBruijnGraph& dbg, int max_reads = -1)
{
	ifstream file(filename);
	if (!file.is_open()) { cerr << "FASTQ 파일 열기 실패: " << filename << "\n"; return 0; }

	string line;
	long long count = 0;
	int line_in_record = 0;

	while (getline(file, line))
	{
		if (line_in_record == 1)
		{
			for (int i = 0; i < (int)line.size(); i++)
				line[i] = DNA_UPPER[(unsigned char)line[i]];
			dbg.add_read(line);
			count++;
			if (max_reads > 0 && count >= max_reads) break;
		}
		line_in_record = (line_in_record + 1) % 4;
	}
	return count;
}

// ──────────────────────────────────────────
// 4. 시뮬레이션
// ──────────────────────────────────────────
void simulate_reads(const string& dna, DeBruijnGraph& dbg, int M, int L)
{
	mt19937 gen(random_device{}());
	uniform_int_distribution<int> dist(0, (int)dna.length() - L);

	string buf(L, ' ');
	for (int i = 0; i < M; i++)
	{
		buf.assign(dna, dist(gen), L);
		dbg.add_read(buf);
	}
}

// ──────────────────────────────────────────
// 5. k-mer 지표
// ──────────────────────────────────────────
struct KmerMetrics { double recall, precision, f1; int orig, asm_, common; };

static HashSet<uint64_t> build_kmer_set(const string& seq, int k)  // uint64_t로 변경
{
	HashSet<uint64_t> s;
	const int n = (int)seq.length();
	if (n < k) return s;
	s.reserve((n - k + 1) * 2);

	// k=32이면 2*32=64비트 → mask가 0이 되므로 전체 비트 사용 (별도 처리)
	const uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1ULL) : ~0ULL;  // uint64_t로 변경
	uint64_t val = 0;
	for (int i = 0; i < k - 1; i++) { val <<= 2; val |= encode_base(seq[i]); }
	for (int i = k - 1; i < n; i++)
	{
		val = ((val << 2) & mask) | encode_base(seq[i]);
		s.insert(val);
	}
	return s;
}

KmerMetrics calculate_kmer_metrics(
	const string& original, const vector<string>& contigs, int k)
{
	HashSet<uint64_t> orig_set = build_kmer_set(original, k);  // uint64_t로 변경

	HashSet<uint64_t> asm_set;
	asm_set.reserve(orig_set.size() * 2);
	const uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1ULL) : ~0ULL;  // uint64_t로 변경

	for (const auto& c : contigs)
	{
		const int n = (int)c.length();
		if (n < k) continue;
		uint64_t val = 0;
		for (int i = 0; i < k - 1; i++) { val <<= 2; val |= encode_base(c[i]); }
		for (int i = k - 1; i < n; i++)
		{
			val = ((val << 2) & mask) | encode_base(c[i]);
			asm_set.insert(val);
		}
	}

	int common = 0;
	for (const auto& km : asm_set)
		if (orig_set.count(km)) common++;

	double recall = orig_set.empty() ? 0.0 : (double)common / orig_set.size() * 100.0;
	double precision = asm_set.empty() ? 0.0 : (double)common / asm_set.size() * 100.0;
	double f1 = (recall + precision > 0)
		? 2.0 * recall * precision / (recall + precision) : 0.0;

	return { recall, precision, f1,
			 (int)orig_set.size(), (int)asm_set.size(), common };
}

// ──────────────────────────────────────────
// 6. N50
// ──────────────────────────────────────────
int calculate_n50(const vector<string>& contigs)
{
	if (contigs.empty()) return 0;
	vector<int> lens;
	lens.reserve(contigs.size());
	int total = 0;
	for (const auto& c : contigs) { lens.push_back((int)c.length()); total += (int)c.length(); }
	sort(lens.rbegin(), lens.rend());
	int cum = 0;
	for (int l : lens) { cum += l; if (cum * 2 >= total) return l; }
	return 0;
}

// ──────────────────────────────────────────
// main
//   실제 FASTQ:              ./program reads.fastq
//   FASTQ + 레퍼런스 지표:   ./program reads.fastq ref.fa
//   시뮬레이션:              ./program
// ──────────────────────────────────────────
int main(int argc, char* argv[])
{
	init_tables();

	auto T0 = chrono::high_resolution_clock::now();

	const int k = 21;
	const int min_weight = 2;
	const int min_contig = 2 * k;

	DeBruijnGraph dbg(k);
	string dna;

	if (argc >= 2)
	{
		cout << "Mode: FASTQ streaming (" << argv[1] << ")\n";
		if (argc >= 3) dna = load_fasta(argv[2]);

		long long n = stream_fastq(argv[1], dbg);
		cout << "Reads processed : " << n << "\n";
	}
	else
	{
		dna = load_fasta("../chr22.fa");
		if (dna.empty()) return 1;

		const int M = 5000000;
		const int L = 100;
		cout << "Mode: simulation (M=" << M << ", L=" << L << ", k=" << k << ")\n";
		simulate_reads(dna, dbg, M, L);
	}

	auto T1 = chrono::high_resolution_clock::now();
	cout << "Read -> Graph   : "
		<< chrono::duration_cast<chrono::milliseconds>(T1 - T0).count() << " ms\n";

	dbg.prune(min_weight);

	vector<string> contigs = dbg.generate_contigs(min_contig);

	auto T2 = chrono::high_resolution_clock::now();
	cout << "Contig Build    : "
		<< chrono::duration_cast<chrono::milliseconds>(T2 - T1).count() << " ms\n";

	int n50 = calculate_n50(contigs);
	size_t longest = 0, total_asm = 0;
	for (const auto& c : contigs)
	{
		longest = max(longest, c.length());
		total_asm += c.length();
	}

	cout << "\n===== Result =====\n";
	cout << "k                : " << k << "\n";
	cout << "Assembly Length  : " << total_asm << "\n";
	cout << "Contig Count     : " << contigs.size() << "\n";
	cout << "Longest Contig   : " << longest << "\n";
	cout << "N50              : " << n50 << "\n";

	if (!dna.empty())
	{
		KmerMetrics metrics = calculate_kmer_metrics(dna, contigs, k);
		auto T3 = chrono::high_resolution_clock::now();
		cout << "Metric          : "
			<< chrono::duration_cast<chrono::milliseconds>(T3 - T2).count() << " ms\n";
		cout << "Original Length  : " << dna.length() << "\n";
		cout << "k-mer Recall     : " << metrics.recall << "%\n";
		cout << "k-mer Precision  : " << metrics.precision << "%\n";
		cout << "k-mer F1         : " << metrics.f1 << "%\n";
	}

	cout << "Total           : "
		<< chrono::duration_cast<chrono::milliseconds>(
			chrono::high_resolution_clock::now() - T0).count() << " ms\n";

	return 0;
}