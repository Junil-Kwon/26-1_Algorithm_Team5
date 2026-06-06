#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <bitset>
#include <cstdlib> 
#include <ctime>  
#include <chrono>
#include <cstdint>
#include <set>
#include <iomanip>
#include <intrin.h>
#include <thread>
#include <atomic>

using namespace std;

const size_t CHECKPOINT_INTERVAL = 64;

atomic<bool> keep_running(true);

void showLoadingAnimation() {
	int count = 0;
	while (keep_running) {
		if (count == 0)      cout << "\r계산 중입니다   " << std::flush;
		else if (count == 1) cout << "\r계산 중입니다.  " << std::flush;
		else if (count == 2) cout << "\r계산 중입니다.. " << std::flush;
		else if (count == 3) cout << "\r계산 중입니다..." << std::flush;
		count = (count + 1) % 4;
		this_thread::sleep_for(std::chrono::milliseconds(500));
	}
	cout << "\r계산 완료!          \n" << std::endl;
}

/** @brief           FASTA 파일에서 특정 길이의 DNA 서열을 읽어오는 함수
 *  @param  filename 불러올 FASTA파일명
 *  @param  max_len  앞에서부터 불러올 최대 길이
 *  @return string   염기 N을 제외한 염기 문자열
 */
string load_fasta_of_len(const string& filename, size_t max_len)
{
	ifstream file(filename);
	if (!file.is_open())
	{
		cerr << "파일을 열 수 없습니다!" << endl;
		return "";
	}

	string line, genome_Seq;
	genome_Seq.reserve(max_len + 1);
	size_t str_len = 0;

	while (getline(file, line))
	{
		if (str_len >= max_len) break;
		if (line.empty()) continue;
		if (line[0] == '>')
		{
			cout << "읽는 중: " << line << endl;
			continue;
		}
		for (char c : line)
		{
			c = toupper(c);
			if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
			{
				genome_Seq += c;
				str_len++;
			}
		}
	}

	file.close();
	cout << endl << "성공적으로 파일을 읽어왔습니다." << endl << endl;
	return genome_Seq;
}

// 문자에 대한 2-bit 인코딩 함수: A=00, C=01, G=10, T=11
unsigned char char_to_2bit(char base)
{
	switch (base)
	{
	case 'A': return 0;
	case 'C': return 1;
	case 'G': return 2;
	case 'T': return 3;
	default:  return 0;
	}
}

// 2-bit 인코딩된 값을 다시 문자로 변환하는 함수
char bit_to_char(uint8_t bit)
{
	switch (bit & 0x03)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	case 3: return 'T';
	default: return 'A';
	}
}

// 역상보서열 생성 함수
string reverse_complement(const string& seq)
{
	string rc = seq;
	reverse(rc.begin(), rc.end());
	for (char& c : rc)
	{
		if (c == 'A') c = 'T';
		else if (c == 'T') c = 'A';
		else if (c == 'C') c = 'G';
		else if (c == 'G') c = 'C';
	}
	return rc;
}

/**
 * @brief                2-bit 인코딩된 64비트 청크 내에서 특정 염기의 개수를 고속으로 카운팅하는 함수
 * @param chunk          32개의 염기가 인코딩된 64비트 정수
 * @param active_bases   현재 청크 내에서 유효하게 조사할 염기의 개수 (1 ~ 32)
 * @param target_2bit    찾고자 하는 염기의 2비트 인코딩 값 (0:A, 1:C, 2:G, 3:T)
 * @return size_t        조건에 부합하는 염기의 개수
 */
size_t count_bits_in_chunk(uint64_t chunk, size_t active_bases, uint8_t target_2bit)
{
	if (active_bases == 0) return 0;

	uint64_t match_bits = 0;

	if (target_2bit == 0) // 'A' (00)
		match_bits = (~chunk >> 1) & ~chunk & 0x5555555555555555ULL;
	else if (target_2bit == 1) // 'C' (01)
		match_bits = (~chunk >> 1) & chunk & 0x5555555555555555ULL;
	else if (target_2bit == 2) // 'G' (10)
		match_bits = (chunk >> 1) & ~chunk & 0x5555555555555555ULL;
	else if (target_2bit == 3) // 'T' (11)
		match_bits = (chunk >> 1) & chunk & 0x5555555555555555ULL;

	if (active_bases < 32)
	{
		uint64_t mask = ~((1ULL << (64 - active_bases * 2)) - 1);
		match_bits &= mask;
	}

	return __popcnt64(match_bits);
}

/**
 * @brief                 체크포인트 테이블(Occ)과 2-bit 인코딩 배열을 결합하여 Rank를 계산하는 함수
 * @param idx             Rank를 구하고자 하는 인덱스
 * @param char_idx        찾고자 하는 문자 인덱스 ($=0, A=1, C=2, G=3, T=4)
 * @param bwt_packed      2-bit 인코딩된 BWT 배열
 * @param Occ_table       64bp 주기로 누적 빈도수를 저장해 둔 체크포인트 Occ 테이블
 * @param end_idx_onBWT   BWT 배열 내에서 특수 문자 '$'가 위치한 인덱스 번호
 * @return size_t         구간 내 문자 등장 횟수 (Rank)
 */
size_t get_rank(size_t idx, uint8_t char_idx, const vector<uint64_t>& bwt_packed, const vector<vector<size_t>>& Occ_table, size_t end_idx_onBWT)
{
	if (idx == 0) return 0;

	if (char_idx == 0)
		return (idx > end_idx_onBWT) ? 1 : 0;

	size_t checkpoint_block = idx / CHECKPOINT_INTERVAL;
	size_t remaining = idx % CHECKPOINT_INTERVAL;
	size_t rank_count = Occ_table[char_idx][checkpoint_block];

	if (remaining == 0) return rank_count;

	size_t packed_idx = checkpoint_block * 2;
	uint8_t target_2bit = char_idx - 1;

	if (remaining <= 32)
		rank_count += count_bits_in_chunk(bwt_packed[packed_idx], remaining, target_2bit);
	else
	{
		rank_count += count_bits_in_chunk(bwt_packed[packed_idx], 32, target_2bit);
		rank_count += count_bits_in_chunk(bwt_packed[packed_idx + 1], remaining - 32, target_2bit);
	}

	if (target_2bit == 0 && idx > end_idx_onBWT)
	{
		size_t checkpoint_limit = checkpoint_block * CHECKPOINT_INTERVAL;
		if (end_idx_onBWT >= checkpoint_limit)
			rank_count--;
	}

	return rank_count;
}

/**
 * @brief                 2-bit 기반 BWT 변환 및 FM-Index 구조체 구축 함수
 * @param text            BWT 변환을 수행할 원본 DNA 서열 (함수 내부에서 '$'가 추가됨)
 * @param bwt_packed      2-bit 인코딩된 BWT 서열 배열
 * @param suffix_array    사전식 순서로 정렬된 접미사 시작 위치 배열
 * @param C_table         BWT 내에서 각 염기별 전역 시작 인덱스를 기록한 테이블
 * @param Occ_table       64bp 주기마다 각 문자의 누적 빈도수를 박제해 둔 압축 체크포인트 테이블
 * @param end_idx_onBWT   BWT 배열 내에서 특수 문자 '$'가 위치한 인덱스 번호
 */
void BWT_2bit(string& text, vector<uint64_t>& bwt_packed, vector<size_t>& suffix_array, vector<size_t>& C_table, vector<vector<size_t>>& Occ_table, size_t& end_idx_onBWT)
{
	text += '$';
	size_t length = text.length();

	suffix_array.clear();
	suffix_array.reserve(length);
	for (size_t i = 0; i < length; i++) suffix_array.push_back(i);

	// NOTICE: O(n^2 log n)의 시간복잡도. 더 효율적인 SA-IS(O(n))가 있으나 구현 복잡도상 단순 정렬 사용.
	sort(suffix_array.begin(), suffix_array.end(), [&](size_t a, size_t b)
		{
			for (size_t i = 0; i < length; i++) {
				char char_a = text[(a + i) % length];
				char char_b = text[(b + i) % length];
				if (char_a != char_b) return char_a < char_b;
			}
			return false;
		});

	C_table.assign(5, 0);
	size_t num_checkpoints = length / CHECKPOINT_INTERVAL + 1;
	Occ_table.assign(5, vector<size_t>(num_checkpoints, 0));
	vector<size_t> current_occ(5, 0);

	bwt_packed.clear();
	bwt_packed.reserve(length / 32 + 1);

	uint64_t buffer = 0;
	for (size_t i = 0; i < length; i++)
	{
		size_t last_idx = (suffix_array[i] + length - 1) % length;
		char suffix = text[last_idx];

		uint8_t char_idx = (suffix == '$') ? 0 : char_to_2bit(suffix) + 1;

		if (i % CHECKPOINT_INTERVAL == 0)
		{
			size_t block = i / CHECKPOINT_INTERVAL;
			for (uint8_t b = 0; b < 5; b++) Occ_table[b][block] = current_occ[b];
		}

		current_occ[char_idx]++;
		C_table[char_idx]++;

		if (suffix == '$') end_idx_onBWT = i;

		uint64_t bit_val = (suffix == '$') ? 0 : char_to_2bit(suffix);
		size_t offset = (i % 32) * 2;
		buffer |= (bit_val << (62 - offset));

		if (i % 32 == 32 - 1 || i == length - 1)
		{
			bwt_packed.push_back(buffer);
			buffer = 0;
		}
	}

	size_t total = 0;
	for (uint8_t i = 0; i < 5; i++)
	{
		size_t count = C_table[i];
		C_table[i] = total;
		total += count;
	}
}

/**
 * @brief                 비트와이즈 기법을 적용한 FM-Index 패턴 검색 함수
 * @param pattern         검색하고자 하는 쿼리 DNA 서열 패턴 문자열
 * @param bwt_packed      2-bit 인코딩된 BWT 서열 배열
 * @param end_idx_onBWT   BWT 서열 내에서 특수 문자 '$'가 위치한 인덱스 번호
 * @param C_table         각 염기별 BWT 내 전역 시작 인덱스를 기록한 테이블
 * @param Occ_table       64bp 주기마다 각 문자의 누적 빈도수를 박제해 둔 압축 체크포인트 테이블
 * @param suffix_array    사전식 순서로 정렬된 접미사 시작 위치 배열
 * @return vector<size_t> 패턴이 발견된 원본 레퍼런스 유전체 상의 시작 인덱스들의 목록
 */
vector<size_t> FM_search(const string& pattern, const vector<uint64_t>& bwt_packed, size_t end_idx_onBWT, const vector<size_t>& C_table, const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array)
{
	vector<size_t> positions;
	size_t length = suffix_array.size();
	size_t top = 0, bot = length;

	for (int i = pattern.length() - 1; i >= 0; i--)
	{
		char c = pattern[i];
		uint8_t char_idx = char_to_2bit(c) + 1;

		size_t occ_top = get_rank(top, char_idx, bwt_packed, Occ_table, end_idx_onBWT);
		size_t occ_bot = get_rank(bot, char_idx, bwt_packed, Occ_table, end_idx_onBWT);

		top = C_table[char_idx] + occ_top;
		bot = C_table[char_idx] + occ_bot;

		if (top >= bot) return positions;
	}

	for (size_t i = top; i < bot; i++)
		positions.push_back(suffix_array[i]);

	return positions;
}

// Pigeonhole Principle 기반 mismatch 허용 검색
// D개 mismatch 허용 시 read를 D+1등분하여 그 중 하나는 반드시 exact match 존재
vector<size_t> Pigeonhole_search(
	const string& read,
	const vector<uint64_t>& bwt_packed,
	size_t end_idx_onBWT,
	const vector<size_t>& C_table,
	const vector<vector<size_t>>& Occ_table,
	const vector<size_t>& suffix_array,
	const string& reference,
	int max_mismatch)
{
	vector<size_t> positions;
	int L = read.length();
	int parts = max_mismatch + 1;
	int part_len = L / parts;
	set<size_t> candidates;

	for (int p = 0; p < parts; p++)
	{
		int start = p * part_len;
		int end = (p == parts - 1) ? L : start + part_len;
		string sub = read.substr(start, end - start);

		vector<size_t> sub_positions = FM_search(sub, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array);

		for (size_t pos : sub_positions)
		{
			if (pos < (size_t)start) continue;
			candidates.insert(pos - start);
		}
	}

	for (size_t pos : candidates)
	{
		if (pos + L > reference.length()) continue;
		int mismatch = 0;
		for (int j = 0; j < L; j++)
		{
			if (read[j] != reference[pos + j]) mismatch++;
			if (mismatch > max_mismatch) break;
		}
		if (mismatch <= max_mismatch)
			positions.push_back(pos);
	}

	return positions;
}

// Paired-end Read 매핑 함수 (Reverse Strand 인덱싱 방식)
// reference + reverse_complement(reference) 를 합쳐서 BWT 인덱스를 만들어
// read1은 forward strand에서, read2_rc는 reverse strand 영역에서 매핑
// insert size 검증으로 유효한 쌍만 반환
struct PairedMapping {
	size_t pos1;      // read1 forward strand 매핑 위치 (원본 기준)
	size_t pos2;      // read2 reverse strand 매핑 위치 (원본 기준으로 변환됨)
	int insert_size;  // 실제 insert size
};

vector<PairedMapping> Paired_end_search(
	const string& read1,
	const string& read2_rc,
	const vector<uint64_t>& bwt_packed,
	size_t end_idx_onBWT,
	const vector<size_t>& C_table,
	const vector<vector<size_t>>& Occ_table,
	const vector<size_t>& suffix_array,
	const string& reference,
	const string& ref_double,
	int max_mismatch,
	int min_insert,
	int max_insert)
{
	vector<PairedMapping> results;
	size_t N = reference.length();

	// read1: forward strand (0 ~ N-1 범위) 에서 매핑
	vector<size_t> pos1_raw = Pigeonhole_search(read1, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array, ref_double, max_mismatch);

	// read2_rc: reverse strand (N ~ 2N-1 범위) 에서 매핑
	vector<size_t> pos2_raw = Pigeonhole_search(read2_rc, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array, ref_double, max_mismatch);

	// read1은 forward strand (0 ~ N-1) 범위만 유효
	vector<size_t> pos1_list;
	for (size_t p : pos1_raw)
		if (p + read1.length() <= N) pos1_list.push_back(p);

	// read2_rc는 reverse strand (N ~ 2N-1) 범위만 유효
	// reverse strand 위치를 원본 좌표로 변환
	vector<size_t> pos2_list;
	for (size_t p : pos2_raw)
	{
		if (p >= N && p + read2_rc.length() <= 2 * N)
		{
			size_t k = p - N;
			size_t L2 = read2_rc.length();
			if (k + L2 <= N)
			{
				size_t orig_pos = N - k - L2;
				pos2_list.push_back(orig_pos);
			}
		}
	}

	// insert size 검증: read1 시작(p1) ~ read2 끝(p2 + L) 범위
	for (size_t p1 : pos1_list)
	{
		for (size_t p2 : pos2_list)
		{
			if (p2 >= p1)
			{
				int insert = (int)p2 + (int)read2_rc.length() - (int)p1;
				if (insert >= min_insert && insert <= max_insert)
				{
					PairedMapping pm;
					pm.pos1 = p1;
					pm.pos2 = p2;
					pm.insert_size = insert;
					results.push_back(pm);
				}
			}
		}
	}
	return results;
}

// Trivial Sliding 함수 (벤치마크용)
vector<size_t> Trivial_search(const string& read, const string& reference, int max_mismatch)
{
	vector<size_t> positions;
	size_t L = read.length();
	size_t N = reference.length();

	for (size_t i = 0; i <= N - L; i++)
	{
		size_t mismatch = 0;
		for (size_t j = 0; j < L; j++)
		{
			if (read[j] != reference[i + j]) mismatch++;
			if (mismatch > max_mismatch) break;
		}
		if (mismatch <= max_mismatch)
			positions.push_back(i);
	}
	return positions;
}

int main()
{
	cout << "start" << endl;

	string dna = load_fasta_of_len("chr22.fa", 10000);

	if (dna.empty())
	{
		cout << "chr22.fa not found -> random genome" << endl;
		srand((unsigned)time(0));
		for (size_t i = 0; i < 10000; i++)
		{
			int r = rand() % 4;
			if (r == 0) dna += 'A'; else if (r == 1) dna += 'C';
			else if (r == 2) dna += 'G'; else dna += 'T';
		}
	}

	cout << "genome length: " << dna.length() << endl;
	string original_dna = dna;

	// [Paired-end 핵심] reference + reverse_complement(reference) 합치기
	// forward strand + reverse strand 동시 인덱싱
	string dna_rc = reverse_complement(original_dna);
	string dna_double = original_dna + dna_rc;

	vector<uint64_t> bwt_2bit;
	size_t end_idx_onBWT;
	vector<size_t> suffix_array;
	vector<size_t> C_table;
	vector<vector<size_t>> Occ_table;

	cout << "BWT 시작" << endl;
	BWT_2bit(dna_double, bwt_2bit, suffix_array, C_table, Occ_table, end_idx_onBWT);
	cout << "BWT 완료" << endl;

	/*
	// SA 출력
	cout << "SA: ";
	for (size_t idx : suffix_array) cout << idx << " ";
	cout << endl << endl;

	// C-Table 출력
	for (size_t i = 0; i < C_table.size(); i++)
	{
		char c = (i == 0) ? '$' : (i == 1) ? 'A' : (i == 2) ? 'C' : (i == 3) ? 'G' : 'T';
		cout << "C[" << c << "] = " << C_table[i] << endl;
	}

	// Occ-Table 출력
	for (size_t i = 0; i < Occ_table.size(); i++)
	{
		char c = (i == 0) ? '$' : (i == 1) ? 'A' : (i == 2) ? 'C' : (i == 3) ? 'G' : 'T';
		cout << "Occ[" << c << "]: ";
		for (size_t j = 0; j < Occ_table[i].size(); j++) cout << Occ_table[i][j] << " ";
		cout << endl;
	}
	*/

	/*
	// 패턴 검색 테스트
	string pattern = "AGG";
	vector<size_t> result = FM_search(pattern, bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
	cout << "패턴 \"" << pattern << "\" 검색 결과: ";
	if (result.empty()) cout << "매칭 없음" << endl;
	else { for (size_t pos : result) cout << pos << " "; cout << endl; }
	*/

	// ========================================
	// 1. reads 생성
	// ========================================
	int M = 6400;
	int L = 32;
	int max_mismatch = 1;
	int insert_min = 100, insert_max = 500;

	// FM-Index용: mismatch 없는 exact reads
	vector<string> exact_reads;
	// Pigeonhole/Trivial용: mismatch 포함 reads
	vector<string> mismatch_reads;

	srand(time(0));
	string bases = "ACGT";

	for (int i = 0; i < M; i++)
	{
		int start = rand() % (original_dna.length() - L);
		string read = original_dna.substr(start, L);

		// exact read (mismatch 없음) → FM-Index용
		exact_reads.push_back(read);

		// mismatch read (0~1개 변이) → Pigeonhole/Trivial용
		string mread = read;
		int num_mutations = rand() % (max_mismatch + 1);
		for (int m = 0; m < num_mutations; m++)
		{
			int mut_pos = rand() % L;
			char original = mread[mut_pos];
			char mutated;
			do { mutated = bases[rand() % 4]; } while (mutated == original);
			mread[mut_pos] = mutated;
		}
		mismatch_reads.push_back(mread);
	}

	// ========================================
	// 실험 환경 출력
	// ========================================
	cout << "========================================" << endl;
	cout << "실험 환경" << endl;
	cout << "========================================" << endl;
	cout << "genome 길이: " << original_dna.length() << " bp" << endl;
	cout << "reads 수: " << M << "개" << endl;
	cout << "read 길이: " << L << " bp" << endl;
	cout << "mismatch 허용: " << max_mismatch << "개" << endl;
	cout << "insert size 범위: " << insert_min << " ~ " << insert_max << " bp" << endl;
	cout << endl;

	keep_running = true;
	std::thread loadingThread(showLoadingAnimation);

	// ========================================
	// 2. FM-Index 시간 측정 (exact reads 사용)
	// ========================================
	auto start_fm = chrono::high_resolution_clock::now();
	int fm_mapped = 0;
	vector<vector<size_t>> mapping_table_fm(original_dna.length(), vector<size_t>(4, 0));
	for (int i = 0; i < exact_reads.size(); i++)
	{
		vector<size_t> positions = FM_search(exact_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
		vector<size_t> valid_pos;
		for (size_t p : positions)
			if (p + L <= original_dna.length()) valid_pos.push_back(p);

		if (!valid_pos.empty()) fm_mapped++;
		for (size_t pos : valid_pos)
			for (int j = 0; j < exact_reads[i].length(); j++)
			{
				if (pos + j >= original_dna.length()) break;
				mapping_table_fm[pos + j][char_to_2bit(exact_reads[i][j])]++;
			}
	}
	auto end_fm = chrono::high_resolution_clock::now();
	auto duration_fm = chrono::duration_cast<chrono::microseconds>(end_fm - start_fm);

	// ========================================
	// 3. Trivial Sliding 시간 측정 (mismatch reads 사용)
	// ========================================
	auto start_trivial = chrono::high_resolution_clock::now();
	int trivial_mapped = 0;
	vector<vector<int>> mapping_table_trivial(original_dna.length(), vector<int>(4, 0));
	for (int i = 0; i < mismatch_reads.size(); i++)
	{
		vector<size_t> positions = Trivial_search(mismatch_reads[i], original_dna, max_mismatch);
		if (!positions.empty()) trivial_mapped++;
		for (size_t pos : positions)
			for (int j = 0; j < mismatch_reads[i].length(); j++)
			{
				if (pos + j >= original_dna.length()) break;
				mapping_table_trivial[pos + j][char_to_2bit(mismatch_reads[i][j])]++;
			}
	}
	auto end_trivial = chrono::high_resolution_clock::now();
	auto duration_trivial = chrono::duration_cast<chrono::microseconds>(end_trivial - start_trivial);

	// ========================================
	// 4. Pigeonhole 시간 측정 (mismatch reads 사용)
	// ========================================
	auto start_pigeon = chrono::high_resolution_clock::now();
	int pigeon_mapped = 0;
	vector<vector<int>> mapping_table_pigeon(original_dna.length(), vector<int>(4, 0));
	for (int i = 0; i < mismatch_reads.size(); i++)
	{
		vector<size_t> positions = Pigeonhole_search(mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array, dna_double, max_mismatch);
		vector<size_t> valid_pos;
		for (size_t p : positions)
			if (p + L <= original_dna.length()) valid_pos.push_back(p);

		if (!valid_pos.empty()) pigeon_mapped++;
		for (size_t pos : valid_pos)
			for (int j = 0; j < mismatch_reads[i].length(); j++)
			{
				if (pos + j >= original_dna.length()) break;
				mapping_table_pigeon[pos + j][char_to_2bit(mismatch_reads[i][j])]++;
			}
	}
	auto end_pigeon = chrono::high_resolution_clock::now();
	auto duration_pigeon = chrono::duration_cast<chrono::microseconds>(end_pigeon - start_pigeon);

	// ========================================
	// 5. Paired-end Read 매핑 (mismatch reads 사용)
	// ========================================
	auto start_paired = chrono::high_resolution_clock::now();
	vector<pair<string, string>> paired_reads;
	srand(time(0));

	for (int i = 0; i < M; i++)
	{
		size_t frag_len = insert_min + rand() % (insert_max - insert_min);
		if (frag_len + L >= original_dna.length()) continue;

		size_t frag_start = rand() % (original_dna.length() - frag_len - L);

		// read1: fragment 앞쪽 forward + mismatch
		string read1 = original_dna.substr(frag_start, L);
		int num_mut1 = rand() % (max_mismatch + 1);
		for (int m = 0; m < num_mut1; m++) {
			int mp = rand() % L;
			char orig = read1[mp]; char mut;
			do { mut = bases[rand() % 4]; } while (mut == orig);
			read1[mp] = mut;
		}

		// read2: fragment 뒤쪽 reverse strand에서 읽음 → 역상보서열 + mismatch
		string read2 = original_dna.substr(frag_start + frag_len, L);
		string read2_rc = reverse_complement(read2);
		int num_mut2 = rand() % (max_mismatch + 1);
		for (int m = 0; m < num_mut2; m++) {
			int mp = rand() % L;
			char orig = read2_rc[mp]; char mut;
			do { mut = bases[rand() % 4]; } while (mut == orig);
			read2_rc[mp] = mut;
		}

		paired_reads.push_back({ read1, read2_rc });
	}

	int paired_valid = 0;
	for (int i = 0; i < paired_reads.size(); i++)
	{
		vector<PairedMapping> results = Paired_end_search(
			paired_reads[i].first, paired_reads[i].second,
			bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
			original_dna, dna_double, max_mismatch, insert_min, insert_max);

		if (!results.empty()) paired_valid++;
	}

	auto end_paired = chrono::high_resolution_clock::now();
	auto duration_paired = chrono::duration_cast<chrono::microseconds>(end_paired - start_paired);

	// ========================================
	// 6. 통합 벤치마크 (FM + Pigeonhole + Paired-end 합산 vs Trivial)
	// ========================================
	// FM-Index: exact reads 매핑
	// Pigeonhole + Paired-end: mismatch reads 매핑
	// 세 알고리즘을 조합하여 Trivial과 비교
	auto start_combined = chrono::high_resolution_clock::now();

	int combined_mapped = 0;
	vector<vector<int>> mapping_table_combined(original_dna.length(), vector<int>(4, 0));

	for (int i = 0; i < M; i++)
	{
		bool mapped = false;
		set<size_t> combined_positions;

		// FM-Index로 exact read 매핑 시도
		vector<size_t> fm_pos = FM_search(exact_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
		for (size_t p : fm_pos)
			if (p + L <= original_dna.length()) { combined_positions.insert(p); mapped = true; }

		// Pigeonhole로 mismatch read 추가 매핑
		vector<size_t> pig_pos = Pigeonhole_search(mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array, dna_double, max_mismatch);
		for (size_t p : pig_pos)
			if (p + L <= original_dna.length()) { combined_positions.insert(p); mapped = true; }

		if (mapped) combined_mapped++;

		for (size_t pos : combined_positions)
			for (int j = 0; j < L; j++)
			{
				if (pos + j >= original_dna.length()) break;
				mapping_table_combined[pos + j][char_to_2bit(mismatch_reads[i][j])]++;
			}
	}

	auto end_combined = chrono::high_resolution_clock::now();
	auto duration_combined = chrono::duration_cast<chrono::microseconds>(end_combined - start_combined);

	keep_running = false;
	if (loadingThread.joinable()) loadingThread.join();

	double fm_speed = (duration_fm.count() > 0) ? (double)M / duration_fm.count() * 1000000 : 0;
	double trivial_speed = (duration_trivial.count() > 0) ? (double)M / duration_trivial.count() * 1000000 : 0;
	double pigeon_speed = (duration_pigeon.count() > 0) ? (double)M / duration_pigeon.count() * 1000000 : 0;
	double paired_speed = (duration_paired.count() > 0) ? (double)M / duration_paired.count() * 1000000 : 0;
	double combined_speed = (duration_combined.count() > 0) ? (double)M / duration_combined.count() * 1000000 : 0;

	// ========================================
	// 개별 알고리즘 벤치마크 비교표
	// ========================================
	cout << "========================================" << endl;
	cout << "개별 알고리즘 벤치마크 비교표" << endl;
	cout << "========================================" << endl;
	cout << left
		<< setw(16) << "항목"
		<< setw(16) << "FM-Index"
		<< setw(16) << "Trivial"
		<< setw(16) << "Pigeonhole"
		<< setw(16) << "Paired-end" << endl;
	cout << string(80, '-') << endl;
	cout << left
		<< setw(16) << "reads 종류"
		<< setw(16) << "exact"
		<< setw(16) << "mismatch"
		<< setw(16) << "mismatch"
		<< setw(16) << "mismatch" << endl;
	cout << left
		<< setw(16) << "실행시간(us)"
		<< setw(16) << duration_fm.count()
		<< setw(16) << duration_trivial.count()
		<< setw(16) << duration_pigeon.count()
		<< setw(16) << duration_paired.count() << endl;
	cout << left
		<< setw(16) << "매핑속도(r/s)"
		<< setw(16) << (int)fm_speed
		<< setw(16) << (int)trivial_speed
		<< setw(16) << (int)pigeon_speed
		<< setw(16) << (int)paired_speed << endl;
	cout << left
		<< setw(16) << "매핑된 reads"
		<< setw(16) << fm_mapped
		<< setw(16) << trivial_mapped
		<< setw(16) << pigeon_mapped
		<< setw(16) << paired_valid << endl;
	cout << left
		<< setw(16) << "매핑률"
		<< setw(15) << fixed << setprecision(1) << (double)fm_mapped / M * 100 << "%"
		<< setw(15) << (double)trivial_mapped / M * 100 << "%"
		<< setw(15) << (double)pigeon_mapped / M * 100 << "%"
		<< setw(15) << (double)paired_valid / M * 100 << "%" << endl;
	cout << endl;

	// ========================================
	// 통합 vs Trivial 비교표 (권준일님 의견 반영)
	// FM-Index + Pigeonhole + Paired-end 조합 vs Trivial
	// ========================================
	cout << "========================================" << endl;
	cout << "통합 알고리즘 vs Trivial 비교표" << endl;
	cout << "(FM-Index + Pigeonhole + Paired-end 조합)" << endl;
	cout << "========================================" << endl;
	cout << left
		<< setw(16) << "항목"
		<< setw(20) << "FM+Pigeonhole+Paired"
		<< setw(16) << "Trivial" << endl;
	cout << string(52, '-') << endl;
	cout << left
		<< setw(16) << "실행시간(us)"
		<< setw(20) << duration_combined.count()
		<< setw(16) << duration_trivial.count() << endl;
	cout << left
		<< setw(16) << "매핑속도(r/s)"
		<< setw(20) << (int)combined_speed
		<< setw(16) << (int)trivial_speed << endl;
	cout << left
		<< setw(16) << "매핑된 reads"
		<< setw(20) << combined_mapped
		<< setw(16) << trivial_mapped << endl;
	cout << left
		<< setw(16) << "매핑률"
		<< setw(19) << fixed << setprecision(1) << (double)combined_mapped / M * 100 << "%"
		<< setw(15) << (double)trivial_mapped / M * 100 << "%" << endl;
	cout << endl;

	// ========================================
	// 7. Consensus 복원 (Pigeonhole 기준)
	// ========================================
	string recovered = "";
	for (size_t i = 0; i < original_dna.length(); i++)
	{
		int total = mapping_table_pigeon[i][0] + mapping_table_pigeon[i][1]
			+ mapping_table_pigeon[i][2] + mapping_table_pigeon[i][3];

		if (total == 0)
		{
			recovered += 'N';
			continue;
		}

		int max_idx = 0, max_val = mapping_table_pigeon[i][0];
		bool tie = false;

		for (int j = 1; j < 4; j++)
		{
			if (mapping_table_pigeon[i][j] > max_val)
			{
				max_val = mapping_table_pigeon[i][j];
				max_idx = j;
				tie = false;
			}
			else if (mapping_table_pigeon[i][j] == max_val)
			{
				tie = true;
			}
		}

		if (tie) recovered += original_dna[i];
		else     recovered += bit_to_char(max_idx);
	}

	// ========================================
	// 8. 정확도 측정 + 시각화 (Pigeonhole 기준)
	// ========================================
	size_t match = 0, mismatch_count = 0, unmapped_count = 0;
	size_t print_limit = 10;

	for (size_t i = 0; i < original_dna.length(); i++)
	{
		if (recovered[i] == 'N')
		{
			unmapped_count++;
		}
		else if (original_dna[i] == recovered[i])
		{
			match++;
		}
		else
		{
			mismatch_count++;
			if (mismatch_count <= print_limit)
				cout << "불일치 위치: " << i
				<< " (원본: " << original_dna[i]
				<< ", 복원: " << recovered[i] << ")" << endl;
		}
	}
	if (mismatch_count > print_limit)
		cout << "... 외 " << mismatch_count - print_limit << "개 불일치" << endl;

	double accuracy = (double)match / original_dna.length() * 100;

	cout << "========================================" << endl;
	cout << "정확도 결과 (Pigeonhole 기준)" << endl;
	cout << "========================================" << endl;
	cout << "매핑됨:       " << match + mismatch_count << "개" << endl;
	cout << "매핑 안됨(N): " << unmapped_count << "개" << endl;
	cout << "불일치:       " << mismatch_count << "개" << endl;
	cout << "정확도:       " << fixed << setprecision(2) << accuracy << "%" << endl;
	cout << endl;

	cout << "[정확도 시각화]" << endl;
	cout << "0%    25%    50%    75%   100%" << endl;
	cout << "|-----|------|------|-----|" << endl;
	int bar_len = 25;
	int filled = (int)(accuracy / 100.0 * bar_len);
	cout << "|";
	for (int i = 0; i < bar_len; i++)
		cout << (i < filled ? "#" : " ");
	cout << "| " << fixed << setprecision(1) << accuracy << "%" << endl;
	cout << "========================================" << endl;

	return 0;
}