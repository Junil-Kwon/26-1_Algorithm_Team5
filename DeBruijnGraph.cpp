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
#include <random>

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

// DNA lookup tables
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

inline uint64_t encode_base(char c)
{
	return BASE_ENC[(unsigned char)c];
}

inline char decode_base(uint64_t x)
{
	return BASE_DEC[x & 3];
}

string decode_kmer(uint64_t value, int length)
{
	string s(length, 'A');
	for (int i = length - 1; i >= 0; i--) {
		s[i] = decode_base(value);
		value >>= 2;
	}
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
	if (max_len > 0) genome.reserve(max_len);

	while (getline(f, line)) {
		if (line.empty() || line[0] == '>') continue;
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

void add_random_mismatch(string& read, int max_mismatch, mt19937& rng)
{
	if (max_mismatch <= 0 || read.empty()) return;

	uniform_int_distribution<int> mutation_count_dist(0, max_mismatch);
	uniform_int_distribution<int> pos_dist(0, (int)read.size() - 1);
	uniform_int_distribution<int> base_dist(0, 3);

	int num_mutations = mutation_count_dist(rng);
	static const char bases[4] = { 'A', 'C', 'G', 'T' };

	for (int m = 0; m < num_mutations; m++) {
		int pos = pos_dist(rng);
		char original = read[pos];
		char mutated;
		do {
			mutated = bases[base_dist(rng)];
		} while (mutated == original);
		read[pos] = mutated;
	}
}

struct PairedRead {
	string read1;
	string read2_rc;
	int insert_size;
};

vector<PairedRead> generate_paired_reads(
	const string& genome,
	int pair_count,
	int read_len,
	int insert_min,
	int insert_max,
	int max_mismatch,
	mt19937& rng)
{
	vector<PairedRead> pairs;
	pairs.reserve(pair_count);

	int n = (int)genome.size();
	if (n <= read_len + insert_min) return pairs;

	uniform_int_distribution<int> frag_len_dist(insert_min, insert_max);

	for (int i = 0; i < pair_count; i++) {
		int frag_len = frag_len_dist(rng);
		if (frag_len + read_len >= n) continue;

		uniform_int_distribution<int> start_dist(0, n - frag_len - read_len);
		int frag_start = start_dist(rng);

		string read1 = genome.substr(frag_start, read_len);
		add_random_mismatch(read1, max_mismatch, rng);

		string read2 = genome.substr(frag_start + frag_len, read_len);
		string read2_rc = reverse_complement(read2);
		add_random_mismatch(read2_rc, max_mismatch, rng);

		PairedRead pr;
		pr.read1 = read1;
		pr.read2_rc = read2_rc;
		pr.insert_size = frag_len + read_len;
		pairs.push_back(pr);
	}

	return pairs;
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
		else node_mask = ~0ULL;
		adj.reserve(1 << 20);
	}

	void add_read(const string& read)
	{
		int n = (int)read.size();
		if (n < K) return;

		uint64_t prefix = 0;
		for (int i = 0; i < K - 1; i++) {
			prefix = ((prefix << 2) & node_mask) | encode_base(read[i]);
		}

		if (adj.find(prefix) == adj.end()) adj[prefix] = vector<Edge>();

		for (int i = K - 1; i < n; i++) {
			uint64_t b = encode_base(read[i]);
			uint64_t suffix = ((prefix << 2) & node_mask) | b;

			vector<Edge>& edges = adj[prefix];
			bool found = false;

			for (size_t j = 0; j < edges.size(); j++) {
				if (edges[j].to == suffix) {
					edges[j].weight++;
					found = true;
					break;
				}
			}

			if (!found) {
				Edge e;
				e.to = suffix;
				e.weight = 1;
				edges.push_back(e);
				if (adj.find(suffix) == adj.end()) adj[suffix] = vector<Edge>();
			}

			prefix = suffix;
		}
	}

	void prune(int min_weight)
	{
		for (auto& kv : adj) {
			vector<Edge>& edges = kv.second;
			edges.erase(
				remove_if(edges.begin(), edges.end(),
					[min_weight](const Edge& e) { return (int)e.weight < min_weight; }),
				edges.end()
			);
		}
	}
};

size_t count_graph_edges(const DeBruijnGraph& dbg)
{
	size_t total = 0;
	for (const auto& kv : dbg.adj) total += kv.second.size();
	return total;
}

vector<string> generate_contigs_unitig(DeBruijnGraph& dbg, int min_len)
{
	HashMap<uint64_t, int> indeg;
	HashMap<uint64_t, int> outdeg;
	indeg.reserve(dbg.adj.size() * 2);
	outdeg.reserve(dbg.adj.size() * 2);

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
	visited.reserve(count_graph_edges(dbg) * 2 + 1);
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

struct KmerIdx {
	int cid;
	int pos;
};

struct KmerHit {
	int cid;
	int pos;
};

HashMap<uint64_t, vector<KmerIdx>> build_kmer_index(const vector<string>& contigs, int klen)
{
	HashMap<uint64_t, vector<KmerIdx>> idx;
	idx.reserve(contigs.size() * 200);

	uint64_t mask = (klen < 32) ? (1ULL << (2 * klen)) - 1ULL : ~0ULL;

	for (int ci = 0; ci < (int)contigs.size(); ci++) {
		const string& c = contigs[ci];
		if ((int)c.size() < klen) continue;

		uint64_t val = 0;
		for (int i = 0; i < klen - 1; i++) {
			val = ((val << 2) & mask) | encode_base(c[i]);
		}

		for (int i = klen - 1; i < (int)c.size(); i++) {
			val = ((val << 2) & mask) | encode_base(c[i]);
			KmerIdx ki;
			ki.cid = ci;
			ki.pos = i - klen + 1;
			idx[val].push_back(ki);
		}
	}

	return idx;
}

HashSet<uint64_t> find_repetitive_kmers(
	const HashMap<uint64_t, vector<KmerIdx>>& idx,
	int max_occurrence)
{
	HashSet<uint64_t> repetitive;
	for (const auto& kv : idx) {
		if ((int)kv.second.size() > max_occurrence) repetitive.insert(kv.first);
	}
	return repetitive;
}

KmerHit map_read_improved(
	const string& read,
	const HashMap<uint64_t, vector<KmerIdx>>& idx,
	const HashSet<uint64_t>& repetitive_kmers,
	int klen,
	const vector<string>& contigs,
	int max_mismatch)
{
	KmerHit fail;
	fail.cid = -1;
	fail.pos = -1;

	if ((int)read.size() < klen) return fail;

	uint64_t mask = (klen < 32) ? (1ULL << (2 * klen)) - 1ULL : ~0ULL;

	// key: (contig id, estimated read start position), value: vote count
	map<pair<int, int>, int> vote;
	int step = max(1, klen / 3);

	for (int offset = 0; offset + klen <= (int)read.size(); offset += step) {
		uint64_t val = 0;
		for (int i = 0; i < klen; i++) {
			val = ((val << 2) & mask) | encode_base(read[offset + i]);
		}

		if (repetitive_kmers.count(val)) continue;

		auto it = idx.find(val);
		if (it == idx.end()) continue;

		for (const KmerIdx& hit : it->second) {
			int estimated_start = hit.pos - offset;
			if (estimated_start < 0) continue;
			vote[make_pair(hit.cid, estimated_start)]++;
		}
	}

	if (vote.empty()) return fail;

	int best_cid = -1;
	int best_pos = -1;
	int best_vote = -1;
	int second_vote = -1;

	for (const auto& kv : vote) {
		int cid = kv.first.first;
		int pos = kv.first.second;
		int v = kv.second;

		if (cid < 0 || cid >= (int)contigs.size()) continue;
		if (pos < 0) continue;
		if (pos + (int)read.size() > (int)contigs[cid].size()) continue;

		if (v > best_vote) {
			second_vote = best_vote;
			best_vote = v;
			best_cid = cid;
			best_pos = pos;
		}
		else if (v > second_vote) {
			second_vote = v;
		}
	}

	if (best_cid == -1) return fail;

	// Require at least two supporting k-mers when the read is long enough.
	if (best_vote < 2 && (int)read.size() >= klen * 2) return fail;

	// Ambiguous repetitive hit.
	if (second_vote == best_vote) return fail;

	int mismatch = 0;
	for (int i = 0; i < (int)read.size(); i++) {
		if (read[i] != contigs[best_cid][best_pos + i]) {
			mismatch++;
			if (mismatch > max_mismatch) return fail;
		}
	}

	KmerHit result;
	result.cid = best_cid;
	result.pos = best_pos;
	return result;
}

struct PELink {
	int from;
	int to;
	int gap_estimate;
	int support;
};

vector<PELink> detect_paired_links_improved(
	const vector<PairedRead>& pairs,
	const vector<string>& contigs,
	const HashMap<uint64_t, vector<KmerIdx>>& idx,
	const HashSet<uint64_t>& repetitive_kmers,
	int klen,
	int read_len,
	int max_mismatch)
{
	map<pair<int, int>, pair<long long, int>> acc;

	for (size_t i = 0; i < pairs.size(); i++) {
		string read2_forward = reverse_complement(pairs[i].read2_rc);

		KmerHit h1 = map_read_improved(
			pairs[i].read1, idx, repetitive_kmers, klen, contigs, max_mismatch);

		KmerHit h2 = map_read_improved(
			read2_forward, idx, repetitive_kmers, klen, contigs, max_mismatch);

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
		PELink lk;
		lk.from = kv.first.first;
		lk.to = kv.first.second;
		lk.gap_estimate = (int)(kv.second.first / kv.second.second);
		lk.support = kv.second.second;
		links.push_back(lk);
	}

	sort(links.begin(), links.end(),
		[](const PELink& a, const PELink& b) { return a.support > b.support; });

	return links;
}

struct ScaffoldResult {
	vector<string> scaffolds;
};

// 여러 scaffold를 N gap으로 구분하여 하나의 최종 복원 서열로 표현
// 주의: scaffold 사이의 실제 순서/염기서열이 확정된 것은 아니며,
//       출력 및 평가 편의를 위해 하나의 문자열로 묶는 용도이다.
string build_final_sequence_from_scaffolds(const vector<string>& scaffolds, int gap_size)
{
	string final_sequence;

	if (scaffolds.empty()) {
		return final_sequence;
	}

	size_t total_len = 0;
	for (const string& s : scaffolds) {
		total_len += s.size();
	}

	final_sequence.reserve(total_len + (scaffolds.size() - 1) * gap_size);

	for (size_t i = 0; i < scaffolds.size(); i++) {
		if (i > 0) {
			final_sequence += string(gap_size, 'N');
		}
		final_sequence += scaffolds[i];
	}

	return final_sequence;
}

ScaffoldResult build_scaffolds_all(
	const vector<string>& contigs,
	const vector<PELink>& links,
	int K,
	int min_support)
{
	ScaffoldResult result;
	int n = (int)contigs.size();
	if (n == 0) return result;

	vector<int> best_out(n, -1);
	vector<int> best_out_gap(n, 0);
	vector<int> best_out_support(n, -1);
	vector<int> best_in(n, -1);
	vector<int> best_in_support(n, -1);

	for (size_t i = 0; i < links.size(); i++) {
		const PELink& lk = links[i];
		if (lk.support < min_support) continue;
		if (lk.from < 0 || lk.from >= n || lk.to < 0 || lk.to >= n) continue;
		if (lk.from == lk.to) continue;

		if (lk.support > best_out_support[lk.from]) {
			best_out[lk.from] = lk.to;
			best_out_gap[lk.from] = lk.gap_estimate;
			best_out_support[lk.from] = lk.support;
		}

		if (lk.support > best_in_support[lk.to]) {
			best_in[lk.to] = lk.from;
			best_in_support[lk.to] = lk.support;
		}
	}

	vector<bool> used(n, false);

	vector<int> starts;
	for (int i = 0; i < n; i++) {
		if (best_in[i] == -1) starts.push_back(i);
	}

	sort(starts.begin(), starts.end(), [&](int a, int b) {
		return contigs[a].size() > contigs[b].size();
		});

	vector<int> all_nodes;
	for (int i = 0; i < n; i++) all_nodes.push_back(i);
	sort(all_nodes.begin(), all_nodes.end(), [&](int a, int b) {
		return contigs[a].size() > contigs[b].size();
		});

	auto build_one_chain = [&](int start) -> string {
		vector<int> chain;
		vector<int> gaps;
		HashSet<int> local_seen;

		int cur = start;
		while (cur != -1) {
			if (cur < 0 || cur >= n) break;
			if (used[cur]) break;
			if (local_seen.count(cur)) break;

			local_seen.insert(cur);
			chain.push_back(cur);
			used[cur] = true;

			int nxt = best_out[cur];
			if (nxt == -1) break;
			if (nxt < 0 || nxt >= n) break;
			if (used[nxt]) break;

			gaps.push_back(max(best_out_gap[cur], -(K - 1)));
			cur = nxt;
		}

		if (chain.empty()) return string("");

		string seq = contigs[chain[0]];
		for (int i = 1; i < (int)chain.size(); i++) {
			int gap = gaps[i - 1];
			if (gap <= 0) {
				int trim_val = min(-gap, K - 1);
				trim_val = min(trim_val, (int)seq.size());
				seq.resize(seq.size() - trim_val);
			}
			else {
				seq += string(gap, 'N');
			}
			seq += contigs[chain[i]];
		}

		return seq;
	};

	// First: chains that have no incoming link.
	for (int start : starts) {
		if (!used[start]) {
			string scaf = build_one_chain(start);
			if (!scaf.empty()) result.scaffolds.push_back(scaf);
		}
	}

	// Second: include every remaining contig as a scaffold, too.
	for (int node : all_nodes) {
		if (!used[node]) {
			string scaf = build_one_chain(node);
			if (!scaf.empty()) result.scaffolds.push_back(scaf);
		}
	}

	sort(result.scaffolds.begin(), result.scaffolds.end(),
		[](const string& a, const string& b) { return a.size() > b.size(); });

	return result;
}

string join_contigs_for_evaluation(const vector<string>& contigs)
{
	string assembly;
	size_t total_len = 0;
	for (size_t i = 0; i < contigs.size(); i++) total_len += contigs[i].size();
	assembly.reserve(total_len + contigs.size());

	for (size_t i = 0; i < contigs.size(); i++) {
		assembly += contigs[i];
		assembly += 'N';
	}
	return assembly;
}

struct KmerMetrics {
	double recall;
	double precision;
	double f1;
};

KmerMetrics evaluate_kmer(const string& original, const string& assembly, int k)
{
	auto make_kmer_set = [&](const string& seq) -> HashSet<uint64_t> {
		HashSet<uint64_t> s;
		int n = (int)seq.size();
		if (n < k) return s;

		uint64_t mask = (k < 32) ? (1ULL << (2 * k)) - 1ULL : ~0ULL;
		uint64_t val = 0;
		int valid_count = 0;
		s.reserve(n);

		for (int i = 0; i < n; i++) {
			char c = seq[i];
			if (c == 'N') {
				val = 0;
				valid_count = 0;
				continue;
			}

			val = ((val << 2) & mask) | encode_base(c);
			valid_count++;

			if (valid_count >= k) s.insert(val);
		}

		return s;
	};

	HashSet<uint64_t> orig = make_kmer_set(original);
	HashSet<uint64_t> asm_set = make_kmer_set(assembly);

	int common = 0;
	for (const auto& x : asm_set) {
		if (orig.count(x)) common++;
	}

	KmerMetrics km;
	km.recall = orig.empty() ? 0.0 : (double)common / orig.size() * 100.0;
	km.precision = asm_set.empty() ? 0.0 : (double)common / asm_set.size() * 100.0;
	km.f1 = (km.recall + km.precision > 0.0)
		? 2.0 * km.recall * km.precision / (km.recall + km.precision)
		: 0.0;

	return km;
}

int calculate_n50(const vector<string>& seqs)
{
	if (seqs.empty()) return 0;

	vector<int> lens;
	int total = 0;

	for (size_t i = 0; i < seqs.size(); i++) {
		lens.push_back((int)seqs[i].size());
		total += lens.back();
	}

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

void print_kmer_bar(const string& label, double value)
{
	const int BAR_LEN = 25;
	int filled = (int)(value / 100.0 * BAR_LEN);
	if (filled < 0) filled = 0;
	if (filled > BAR_LEN) filled = BAR_LEN;

	cout << "  " << left << setw(20) << label << " |";
	for (int i = 0; i < BAR_LEN; i++) {
		cout << (i < filled ? "#" : " ");
	}
	cout << "| " << fixed << setprecision(1) << value << " %\n";
}

int main(int argc, char* argv[])
{
	init_tables();

	mt19937 rng((unsigned)chrono::steady_clock::now().time_since_epoch().count());
	cout << fixed << setprecision(4);

	const int K = 21;
	const int PAIR_COUNT = 75000;
	const int READ_LEN = 32;
	const int INSERT_MIN = 150;
	const int INSERT_MAX = 350;
	const int MAX_MISMATCH = 1;
	const int MIN_WEIGHT = 3;
	const int MIN_CONTIG = 100;
	const int MIN_PE_SUPPORT = 2;
	const int MAX_REPEAT_KMER_OCC = 20;
	const int FINAL_GAP_SIZE = 1;
	const size_t GENOME_LEN = 100000;

	string genome_path = (argc >= 2) ? argv[1] : "../chr22.fa";

	auto total_start = chrono::high_resolution_clock::now();

	string genome = load_fasta(genome_path, GENOME_LEN);
	if (genome.empty()) {
		cerr << "genome 로드 실패: " << genome_path << "\n";
		return 1;
	}

	cout << "========================================\n";
	cout << "실험 환경\n";
	cout << "========================================\n";
	cout << "genome 길이     : " << genome.size() << " bp\n";
	cout << "reads 수        : " << PAIR_COUNT * 2 << "개 (paired-end " << PAIR_COUNT << "쌍)\n";
	cout << "read 길이       : " << READ_LEN << " bp\n";
	cout << "k-mer 길이      : " << K << "\n";
	cout << "mismatch 허용   : " << MAX_MISMATCH << "개\n";
	cout << "insert size 범위: " << INSERT_MIN << " ~ " << INSERT_MAX << " bp\n";
	cout << "\n";

	// 1. paired-end read 생성
	vector<PairedRead> paired_reads = generate_paired_reads(
		genome, PAIR_COUNT, READ_LEN, INSERT_MIN, INSERT_MAX, MAX_MISMATCH, rng);

	// 2. De Bruijn Graph 구성 및 pruning
	DeBruijnGraph dbg(K);
	for (size_t i = 0; i < paired_reads.size(); i++) {
		dbg.add_read(paired_reads[i].read1);
		string read2_forward = reverse_complement(paired_reads[i].read2_rc);
		dbg.add_read(read2_forward);
	}
	dbg.prune(MIN_WEIGHT);

	// 3. contig 생성
	vector<string> contigs = generate_contigs_unitig(dbg, MIN_CONTIG);

	// 4. contig k-mer index 생성 및 paired-end link 탐지
	HashMap<uint64_t, vector<KmerIdx>> kmer_idx = build_kmer_index(contigs, K);
	HashSet<uint64_t> repetitive_kmers = find_repetitive_kmers(kmer_idx, MAX_REPEAT_KMER_OCC);
	vector<PELink> links = detect_paired_links_improved(
		paired_reads, contigs, kmer_idx, repetitive_kmers, K, READ_LEN, MAX_MISMATCH);

	// 5. scaffold 생성 후 하나의 최종 복원 서열로 표현
	ScaffoldResult final_result = build_scaffolds_all(contigs, links, K, MIN_PE_SUPPORT);
	string final_sequence = build_final_sequence_from_scaffolds(final_result.scaffolds, FINAL_GAP_SIZE);

	// 6. 정확도 평가
	KmerMetrics final_km = evaluate_kmer(genome, final_sequence, K);

	auto total_end = chrono::high_resolution_clock::now();
	long long runtime_total = elapsed_ms(total_start, total_end);

	size_t total_contig = 0;
	size_t longest_contig = 0;
	for (size_t i = 0; i < contigs.size(); i++) {
		total_contig += contigs[i].size();
		if (contigs[i].size() > longest_contig) longest_contig = contigs[i].size();
	}

	int contig_n50 = calculate_n50(contigs);
	double memory_mb = get_memory_usage_mb();

	cout << "\n========================================\n";
	cout << "De Bruijn Graph Assembly Result\n";
	cout << "========================================\n";

	cout << "\n[자원 사용량]\n";
	cout << "  메모리 사용량        : " << memory_mb << " MB\n";

	cout << "\n[Contig 조립 결과]\n";
	cout << "  Contig 개수          : " << contigs.size() << "\n";
	cout << "  전체 Contig 길이     : " << total_contig << " bp\n";
	cout << "  최장 Contig 길이     : " << longest_contig << " bp\n";
	cout << "  N50                  : " << contig_n50 << " bp\n";

	cout << "\n[최종 복원 서열]\n";
	cout << "  최종 서열 길이       : " << final_sequence.size() << " bp\n";

	cout << "\n[정확도]\n";
	cout << "  " << string(48, '-') << "\n";
	cout << "  " << left << setw(20) << "" << "  0%    25%   50%   75%  100%\n";
	print_kmer_bar("k-mer 재현율", final_km.recall);
	print_kmer_bar("k-mer 정밀도", final_km.precision);
	print_kmer_bar("k-mer F1 점수", final_km.f1);
	cout << "  " << string(48, '-') << "\n";

	cout << "\n[실행시간]\n";
	cout << "  전체 실행시간        : " << runtime_total << " ms\n";

	cout << "========================================\n\n";

	return 0;
}
