#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cctype>

using namespace std;

// FASTA 파일 읽기
string load_fasta(const string& filename)
{
	ifstream file(filename);

	if (!file.is_open())
	{
		cerr << "파일 열기 실패" << endl;
		return "";
	}

	string line;
	string genome;

	while (getline(file, line))
	{
		if (line.empty()) continue;

		if (line[0] == '>')
		{
			continue;
		}

		for (char c : line)
		{
			c = toupper(c);

			// N은 저장하지 않음
			if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
			{
				genome += c;
			}
		}
	}

	file.close();

	return genome;
}

// DNA 문자 -> 2bit
inline uint64_t encode_base(char c)
{
	switch (c)
	{
	case 'A': return 0ULL;
	case 'C': return 1ULL;
	case 'G': return 2ULL;
	case 'T': return 3ULL;
	default: return 0ULL;
	}
}

// 2bit -> DNA 문자
inline char decode_base(uint64_t x)
{
	switch (x & 3ULL)
	{
	case 0ULL: return 'A';
	case 1ULL: return 'C';
	case 2ULL: return 'G';
	case 3ULL: return 'T';
	default: return 'A';
	}
}

// encoded k-mer -> 문자열
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

// 문자열 k-mer -> uint64_t
uint64_t encode_kmer(const string& s)
{
	uint64_t value = 0;

	for (char c : s)
	{
		value <<= 2;
		value |= encode_base(c);
	}

	return value;
}

class DeBruijnGraph
{
private:
	// prefix -> suffix -> weight
	unordered_map<uint64_t, unordered_map<uint64_t, int>> graph;

	// degree는 중복 edge 개수가 아니라 unique edge 기준
	unordered_map<uint64_t, int> indegree;
	unordered_map<uint64_t, int> outdegree;

	int k;
	uint64_t mask;

public:
	DeBruijnGraph(int kmer_size)
	{
		k = kmer_size;

		// uint64_t는 최대 32-mer 정도까지 안전
		if (k > 32)
		{
			cerr << "k가 너무 큽니다. uint64_t 2bit encoding에서는 k <= 32 권장" << endl;
			exit(1);
		}

		mask = (1ULL << (2 * (k - 1))) - 1;
	}

	// read 하나를 rolling 2bit 방식으로 graph에 추가
	void add_read(const string& read)
	{
		if (read.length() < k) return;

		uint64_t prefix = 0;

		// 첫 prefix, 길이 k-1 생성
		for (int i = 0; i < k - 1; i++)
		{
			prefix <<= 2;
			prefix |= encode_base(read[i]);
		}

		// rolling k-mer 생성
		for (int i = k - 1; i < read.length(); i++)
		{
			uint64_t suffix =
				((prefix << 2) & mask) | encode_base(read[i]);

			// 처음 생긴 edge인지 확인
			bool new_edge =
				(graph[prefix].find(suffix) == graph[prefix].end());

			// edge weight 증가
			graph[prefix][suffix]++;

			// degree는 unique edge일 때만 증가
			if (new_edge)
			{
				outdegree[prefix]++;
				indegree[suffix]++;
			}

			// suffix 노드도 graph에 등록
			if (graph.find(suffix) == graph.end())
			{
				graph[suffix] = unordered_map<uint64_t, int>();
			}

			prefix = suffix;
		}
	}

	void print_graph()
	{
		cout << "===== De Bruijn Graph =====" << endl;

		for (auto& node : graph)
		{
			cout << decode_kmer(node.first, k - 1) << " -> ";

			for (auto& edge : node.second)
			{
				cout << decode_kmer(edge.first, k - 1)
					<< "(" << edge.second << ") ";
			}

			cout << endl;
		}

		cout << endl;
	}

	// 가장 weight가 큰 edge 선택
	uint64_t get_best_next(uint64_t node)
	{
		uint64_t best_next = 0;
		int best_weight = -1;

		for (auto& edge : graph[node])
		{
			if (edge.second > best_weight)
			{
				best_weight = edge.second;
				best_next = edge.first;
			}
		}

		return best_next;
	}

	vector<string> generate_contigs()
	{
		vector<string> contigs;

		for (auto& node : graph)
		{
			uint64_t start = node.first;

			// branch 시작점 또는 끝점에서 contig 시작
			if (!(indegree[start] == 1 && outdegree[start] == 1))
			{
				for (auto& edge : graph[start])
				{
					uint64_t current = edge.first;

					string contig = decode_kmer(start, k - 1);
					contig += decode_base(current);

					// 직선 구간이면 계속 진행
					while (indegree[current] == 1 && outdegree[current] == 1)
					{
						uint64_t next_node = get_best_next(current);

						current = next_node;

						contig += decode_base(current);
					}

					contigs.push_back(contig);
				}
			}
		}

		// graph가 완전한 cycle이면 위 조건에서 contig가 안 생길 수 있음
		if (contigs.empty() && !graph.empty())
		{
			uint64_t start = graph.begin()->first;
			uint64_t current = start;

			string contig = decode_kmer(start, k - 1);

			do
			{
				uint64_t next_node = get_best_next(current);
				current = next_node;
				contig += decode_base(current);
			} while (current != start);

			contigs.push_back(contig);
		}

		return contigs;
	}

	int get_overlap(const string& a, const string& b)
	{
		int max_overlap = 0;
		int min_len = min(a.length(), b.length());

		for (int len = 1; len < min_len; len++)
		{
			if (a.substr(a.length() - len) == b.substr(0, len))
			{
				max_overlap = len;
			}
		}

		return max_overlap;
	}

	string merge_contigs(vector<string>& contigs)
	{
		if (contigs.empty()) return "";

		vector<bool> used(contigs.size(), false);

		string merged = contigs[0];
		used[0] = true;

		while (true)
		{
			int best_idx = -1;
			int best_overlap = 0;

			for (int i = 0; i < contigs.size(); i++)
			{
				if (used[i]) continue;

				int overlap = get_overlap(merged, contigs[i]);

				if (overlap > best_overlap)
				{
					best_overlap = overlap;
					best_idx = i;
				}
			}

			if (best_idx == -1) break;

			merged += contigs[best_idx].substr(best_overlap);
			used[best_idx] = true;
		}

		return merged;
	}
};

int main()
{
	/*
	// FASTA 파일 읽기
	string dna = load_fasta("dna.txt");

	if (dna.empty())
	{
		return 1;
	}
	*/

	// 테스트용 DNA
	string dna =
		"AGCTACCAGGTGAGCTTTTTTCCTAGCTTACGAAACGCCCAAATTTTTTAAAAAACCCGAGCGCGAGAGACTACGATGCA";

	cout << "DNA Length: " << dna.length() << endl << endl;

	int read_length = 30;

	// k-mer 길이
	int k = read_length * 2 / 3;

	cout << "Read Length: " << read_length << endl;
	cout << "k-mer Length: " << k << endl << endl;

	vector<string> reads;

	// 테스트용 read 생성
	for (int i = 0; i <= dna.length() - read_length; i++)
	{
		reads.push_back(dna.substr(i, read_length));
	}

	DeBruijnGraph dbg(k);

	for (string& r : reads)
	{
		dbg.add_read(r);
	}

	dbg.print_graph();

	vector<string> contigs = dbg.generate_contigs();

	cout << "===== Contigs =====" << endl;

	for (string& c : contigs)
	{
		cout << c << endl;
	}

	cout << endl;

	string assembled = dbg.merge_contigs(contigs);

	cout << "===== Final Assembly =====" << endl;
	cout << assembled << endl << endl;

	cout << "===== Original DNA =====" << endl;
	cout << dna << endl;

	return 0;
}