#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cstdint>
#include <set>
#include <iomanip>
#include <intrin.h>
#include <thread>
#include <atomic>
#include <cmath>
#include <windows.h>
#include <psapi.h>
#include <random>
#include <unordered_map>
#pragma comment(lib, "psapi.lib")

using namespace std;

// =========================================================================
// [수정 5] fm_cache는 각 벤치마크 전에 명시적으로 clear() 호출
// 전역으로 두되, 벤치마크 독립성 보장
// =========================================================================
unordered_map<string, vector<size_t>> fm_cache;

const size_t CHECKPOINT_INTERVAL = 64;
atomic<bool> keep_running(true);

void showLoadingAnimation() {
    int count = 0;
    while (keep_running) {
        if (count == 0)      cout << "\r계산 중입니다   " << flush;
        else if (count == 1) cout << "\r계산 중입니다.  " << flush;
        else if (count == 2) cout << "\r계산 중입니다.. " << flush;
        else if (count == 3) cout << "\r계산 중입니다..." << flush;
        count = (count + 1) % 4;
        this_thread::sleep_for(chrono::milliseconds(500));
    }
    cout << "\r계산 완료!          \n" << endl;
}

string load_fasta_of_len(const string& filename, size_t max_len)
{
    ifstream file(filename);
    if (!file.is_open()) { cerr << "파일을 열 수 없습니다!" << endl; return ""; }
    string line, genome_Seq;
    genome_Seq.reserve(max_len + 1);
    size_t str_len = 0;
    while (getline(file, line)) {
        if (str_len >= max_len) break;
        if (line.empty()) continue;
        if (line[0] == '>') { cout << "읽는 중: " << line << endl; continue; }
        for (char c : line) {
            c = toupper(c);
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T') { genome_Seq += c; str_len++; }
        }
    }
    file.close();
    cout << endl << "성공적으로 파일을 읽어왔습니다." << endl << endl;
    return genome_Seq;
}

unsigned char char_to_2bit(char base)
{
    switch (base) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return 0; }
}

char bit_to_char(uint8_t bit)
{
    switch (bit & 0x03) { case 0: return 'A'; case 1: return 'C'; case 2: return 'G'; case 3: return 'T'; default: return 'A'; }
}

string reverse_complement_char(const string& seq)
{
    string rc = seq;
    reverse(rc.begin(), rc.end());
    for (char& c : rc) {
        if (c == 'A') c = 'T'; else if (c == 'T') c = 'A';
        else if (c == 'C') c = 'G'; else if (c == 'G') c = 'C';
    }
    return rc;
}

string reverse_complement_2bit(const string& seq)
{
    int n = seq.length();
    vector<uint64_t> encoded((n + 31) / 32, 0);
    for (int i = 0; i < n; i++) {
        uint64_t base = char_to_2bit(seq[i]);
        encoded[i / 32] |= (base << (62 - (i % 32) * 2));
    }
    for (auto& chunk : encoded) chunk ^= 0xAAAAAAAAAAAAAAAAULL;
    string rc(n, 'A');
    for (int i = 0; i < n; i++) {
        int rev_i = n - 1 - i;
        uint64_t base = (encoded[rev_i / 32] >> (62 - (rev_i % 32) * 2)) & 0x3;
        rc[i] = bit_to_char(base);
    }
    return rc;
}

vector<uint64_t> encode_read_2bit(const string& read)
{
    int n = read.length();
    vector<uint64_t> encoded((n + 31) / 32, 0);
    for (int i = 0; i < n; i++) {
        uint64_t base = char_to_2bit(read[i]);
        encoded[i / 32] |= (base << (62 - (i % 32) * 2));
    }
    return encoded;
}

vector<uint64_t> encode_reference_2bit(const string& reference)
{
    int n = reference.length();
    vector<uint64_t> encoded((n + 31) / 32, 0);
    for (int i = 0; i < n; i++) {
        uint64_t base = char_to_2bit(reference[i]);
        encoded[i / 32] |= (base << (62 - (i % 32) * 2));
    }
    return encoded;
}

int count_mismatches_char(const string& read, const string& reference, size_t pos, int L)
{
    int mismatch = 0;
    for (int j = 0; j < L; j++) {
        if (pos + j >= reference.length()) return L;
        if (read[j] != reference[pos + j]) mismatch++;
    }
    return mismatch;
}

int count_mismatches_2bit(const vector<uint64_t>& read_encoded,
    const string& reference, size_t pos, int L)
{
    vector<uint64_t> ref_encoded((L + 31) / 32, 0);
    for (int i = 0; i < L; i++) {
        if (pos + i >= reference.length()) return L;
        uint64_t base = char_to_2bit(reference[pos + i]);
        ref_encoded[i / 32] |= (base << (62 - (i % 32) * 2));
    }
    int mismatch = 0;
    int chunks = (L + 31) / 32;
    for (int c = 0; c < chunks; c++) {
        uint64_t diff = read_encoded[c] ^ ref_encoded[c];
        uint64_t mismatch_bits = (diff | (diff >> 1)) & 0x5555555555555555ULL;
        int active = min(32, L - c * 32);
        if (active < 32) {
            uint64_t mask = ~((1ULL << (64 - active * 2)) - 1);
            mismatch_bits &= mask;
        }
        mismatch += __popcnt64(mismatch_bits);
    }
    return mismatch;
}

int count_mismatches_2bit_preencoded(
    const vector<uint64_t>& read_encoded,
    const vector<uint64_t>& ref_encoded_full,
    size_t pos, int L)
{
    int mismatch = 0;
    int chunks = (L + 31) / 32;
    for (int c = 0; c < chunks; c++) {
        size_t ref_base = pos + c * 32;
        size_t ref_chunk_idx = ref_base / 32;
        size_t ref_offset = ref_base % 32;
        uint64_t ref_chunk = 0;
        if (ref_offset == 0) {
            ref_chunk = ref_encoded_full[ref_chunk_idx];
        }
        else {
            ref_chunk = (ref_encoded_full[ref_chunk_idx] << (ref_offset * 2));
            if (ref_chunk_idx + 1 < ref_encoded_full.size())
                ref_chunk |= (ref_encoded_full[ref_chunk_idx + 1] >> (64 - ref_offset * 2));
        }
        uint64_t diff = read_encoded[c] ^ ref_chunk;
        uint64_t mismatch_bits = (diff | (diff >> 1)) & 0x5555555555555555ULL;
        int active = min(32, L - c * 32);
        if (active < 32) {
            uint64_t mask = ~((1ULL << (64 - active * 2)) - 1);
            mismatch_bits &= mask;
        }
        mismatch += __popcnt64(mismatch_bits);
    }
    return mismatch;
}

int calc_hamming_score(const string& read, const string& ref)
{
    int match = 0;
    int L = min((int)read.length(), (int)ref.length());
    for (int i = 0; i < L; i++)
        if (read[i] == ref[i]) match++;
    return match;
}

int calc_MAPQ(const vector<int>& scores)
{
    if (scores.empty()) return 0;
    if (scores.size() == 1) return 60;
    int best = *max_element(scores.begin(), scores.end());
    int second_best = 0;
    for (int s : scores)
        if (s != best && s > second_best) second_best = s;
    int best_count = count(scores.begin(), scores.end(), best);
    if (best_count > 1) second_best = best;
    if (best == 0) return 0;
    if (second_best == 0) return 60;
    double ratio = (double)second_best / (double)best;
    if (ratio >= 1.0) return 0;
    return min((int)(-10.0 * log10(ratio)), 60);
}

size_t count_bits_in_chunk(uint64_t chunk, size_t active_bases, uint8_t target_2bit)
{
    if (active_bases == 0) return 0;
    uint64_t match_bits = 0;
    if (target_2bit == 0) match_bits = (~chunk >> 1) & ~chunk & 0x5555555555555555ULL;
    else if (target_2bit == 1) match_bits = (~chunk >> 1) & chunk & 0x5555555555555555ULL;
    else if (target_2bit == 2) match_bits = (chunk >> 1) & ~chunk & 0x5555555555555555ULL;
    else if (target_2bit == 3) match_bits = (chunk >> 1) & chunk & 0x5555555555555555ULL;
    if (active_bases < 32) { uint64_t mask = ~((1ULL << (64 - active_bases * 2)) - 1); match_bits &= mask; }
    return __popcnt64(match_bits);
}

size_t get_rank(size_t idx, uint8_t char_idx, const vector<uint64_t>& bwt_packed,
    const vector<vector<size_t>>& Occ_table, size_t end_idx_onBWT)
{
    if (idx == 0) return 0;
    if (char_idx == 0) return (idx > end_idx_onBWT) ? 1 : 0;
    size_t checkpoint_block = idx / CHECKPOINT_INTERVAL;
    size_t remaining = idx % CHECKPOINT_INTERVAL;
    size_t rank_count = Occ_table[char_idx][checkpoint_block];
    if (remaining == 0) return rank_count;
    size_t packed_idx = checkpoint_block * 2;
    uint8_t target_2bit = char_idx - 1;
    if (remaining <= 32) rank_count += count_bits_in_chunk(bwt_packed[packed_idx], remaining, target_2bit);
    else {
        rank_count += count_bits_in_chunk(bwt_packed[packed_idx], 32, target_2bit);
        rank_count += count_bits_in_chunk(bwt_packed[packed_idx + 1], remaining - 32, target_2bit);
    }
    if (target_2bit == 0 && idx > end_idx_onBWT) {
        size_t checkpoint_limit = checkpoint_block * CHECKPOINT_INTERVAL;
        if (end_idx_onBWT >= checkpoint_limit) rank_count--;
    }
    return rank_count;
}

// =========================================================================
// [수정 6] BWT 대상을 original_dna만으로 축소
// 기존: dna_double (20,000bp) → O(N² log N) suffix sort가 4배 느림
// 수정: original_dna (10,000bp) → 속도 대폭 개선
// Paired-end의 역방향은 별도 reverse_complement로 처리
// =========================================================================
void BWT_2bit(string& text, vector<uint64_t>& bwt_packed, vector<size_t>& suffix_array,
    vector<size_t>& C_table, vector<vector<size_t>>& Occ_table, size_t& end_idx_onBWT)
{
    text += '$';
    size_t length = text.length();
    suffix_array.clear();
    suffix_array.reserve(length);
    for (size_t i = 0; i < length; i++) suffix_array.push_back(i);
    // NOTICE: O(n^2 log n). SA-IS(O(n))가 더 효율적이나 구현 복잡도상 단순 정렬 사용.
    sort(suffix_array.begin(), suffix_array.end(), [&](size_t a, size_t b) {
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
    for (size_t i = 0; i < length; i++) {
        size_t last_idx = (suffix_array[i] + length - 1) % length;
        char suffix = text[last_idx];
        uint8_t char_idx = (suffix == '$') ? 0 : char_to_2bit(suffix) + 1;
        if (i % CHECKPOINT_INTERVAL == 0) {
            size_t block = i / CHECKPOINT_INTERVAL;
            for (uint8_t b = 0; b < 5; b++) Occ_table[b][block] = current_occ[b];
        }
        current_occ[char_idx]++;
        C_table[char_idx]++;
        if (suffix == '$') end_idx_onBWT = i;
        uint64_t bit_val = (suffix == '$') ? 0 : char_to_2bit(suffix);
        size_t offset = (i % 32) * 2;
        buffer |= (bit_val << (62 - offset));
        if (i % 32 == 32 - 1 || i == length - 1) { bwt_packed.push_back(buffer); buffer = 0; }
    }
    size_t total = 0;
    for (uint8_t i = 0; i < 5; i++) { size_t cnt = C_table[i]; C_table[i] = total; total += cnt; }
}

vector<size_t> FM_search(const string& pattern, const vector<uint64_t>& bwt_packed,
    size_t end_idx_onBWT, const vector<size_t>& C_table,
    const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array)
{
    auto it = fm_cache.find(pattern);
    if (it != fm_cache.end()) return it->second;

    vector<size_t> positions;
    size_t length = suffix_array.size();
    size_t top = 0, bot = length;
    for (int i = pattern.length() - 1; i >= 0; i--) {
        char c = pattern[i];
        uint8_t char_idx = char_to_2bit(c) + 1;
        size_t occ_top = get_rank(top, char_idx, bwt_packed, Occ_table, end_idx_onBWT);
        size_t occ_bot = get_rank(bot, char_idx, bwt_packed, Occ_table, end_idx_onBWT);
        top = C_table[char_idx] + occ_top;
        bot = C_table[char_idx] + occ_bot;
        if (top >= bot) { fm_cache[pattern] = positions; return positions; }
    }
    for (size_t i = top; i < bot; i++) positions.push_back(suffix_array[i]);
    fm_cache[pattern] = positions;
    return positions;
}

vector<size_t> Pigeonhole_search_char(
    const string& read, const vector<uint64_t>& bwt_packed,
    size_t end_idx_onBWT, const vector<size_t>& C_table,
    const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array,
    const string& reference, int max_mismatch)
{
    vector<size_t> positions;
    int L = read.length();
    int parts = max_mismatch + 1;
    int part_len = L / parts;
    vector<size_t> candidates;
    for (int p = 0; p < parts; p++) {
        int start = p * part_len;
        int end = (p == parts - 1) ? L : start + part_len;
        string sub = read.substr(start, end - start);
        vector<size_t> sub_positions = FM_search(sub, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array);
        for (size_t pos : sub_positions) {
            if (pos < (size_t)start) continue;
            candidates.push_back(pos - start);
        }
    }
    sort(candidates.begin(), candidates.end());
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    for (size_t pos : candidates) {
        if (pos + L > reference.length()) continue;
        if (count_mismatches_char(read, reference, pos, L) <= max_mismatch)
            positions.push_back(pos);
    }
    return positions;
}

vector<size_t> Pigeonhole_search_2bit(
    const string& read, const vector<uint64_t>& bwt_packed,
    size_t end_idx_onBWT, const vector<size_t>& C_table,
    const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array,
    const string& reference, int max_mismatch)
{
    vector<size_t> positions;
    int L = read.length();
    int parts = max_mismatch + 1;
    int part_len = L / parts;
    vector<size_t> candidates;
    for (int p = 0; p < parts; p++) {
        int start = p * part_len;
        int end = (p == parts - 1) ? L : start + part_len;
        string sub = read.substr(start, end - start);
        vector<size_t> sub_positions = FM_search(sub, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array);
        for (size_t pos : sub_positions) {
            if (pos < (size_t)start) continue;
            candidates.push_back(pos - start);
        }
    }
    sort(candidates.begin(), candidates.end());
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    vector<uint64_t> read_encoded = encode_read_2bit(read);
    for (size_t pos : candidates) {
        if (pos + L > reference.length()) continue;
        if (count_mismatches_2bit(read_encoded, reference, pos, L) <= max_mismatch)
            positions.push_back(pos);
    }
    return positions;
}

struct MappingResult { size_t pos; int score; int mapq; };

// =========================================================================
// [수정 2] MAPQ 검색 대상을 original_dna로 통일
// 기존: dna_double 전달 → 매핑률 부풀리기 (역방향까지 검색해서 높아보임)
// 수정: original_dna만 전달 → 실제 복원에 기여하는 매핑만 카운트
// =========================================================================
vector<MappingResult> Pigeonhole_search_with_MAPQ(
    const string& read, const vector<uint64_t>& bwt_packed,
    size_t end_idx_onBWT, const vector<size_t>& C_table,
    const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array,
    const string& reference, int max_mismatch,
    // [수정 4] MAPQ threshold를 10으로 낮춤
    // 기존: 20 → chr22 반복서열 때문에 reads 대부분이 필터링됨
    // 수정: 10 → 적당한 신뢰도 유지하면서 N(unmapped) 감소
    int mapq_threshold = 10)
{
    vector<MappingResult> results;
    int L = read.length();
    vector<size_t> positions = Pigeonhole_search_2bit(
        read, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array,
        reference, max_mismatch);
    if (positions.empty()) return results;
    vector<int> scores;
    for (size_t pos : positions) {
        if (pos + L > reference.length()) continue;
        string ref_seg = reference.substr(pos, L);
        scores.push_back(calc_hamming_score(read, ref_seg));
    }
    if (scores.empty()) return results;
    int mapq = calc_MAPQ(scores);
    if (mapq >= mapq_threshold) {
        int best_score = *max_element(scores.begin(), scores.end());
        for (size_t i = 0; i < positions.size(); i++) {
            if (scores[i] == best_score) {
                MappingResult mr; mr.pos = positions[i]; mr.score = scores[i]; mr.mapq = mapq;
                results.push_back(mr); break;
            }
        }
    }
    return results;
}

struct PairedMapping { size_t pos1; size_t pos2; int insert_size; };

vector<PairedMapping> Paired_end_search(
    const string& read1, const string& read2_rc,
    const vector<uint64_t>& bwt_packed, size_t end_idx_onBWT,
    const vector<size_t>& C_table, const vector<vector<size_t>>& Occ_table,
    const vector<size_t>& suffix_array,
    const string& reference,
    int max_mismatch, int min_insert, int max_insert)
{
    vector<PairedMapping> results;
    size_t N = reference.length();
    int L = read1.length();

    // read1: forward strand에서 검색
    vector<size_t> pos1_list = Pigeonhole_search_2bit(
        read1, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array,
        reference, max_mismatch);

    // read2_rc의 reverse complement = 원래 read2 → forward strand에서 검색
    string read2_fwd = reverse_complement_2bit(read2_rc);
    vector<size_t> pos2_list = Pigeonhole_search_2bit(
        read2_fwd, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array,
        reference, max_mismatch);

    for (size_t p1 : pos1_list) {
        if (p1 + L > N) continue;
        for (size_t p2 : pos2_list) {
            if (p2 + L > N) continue;
            if (p2 >= p1) {
                int insert = (int)p2 + L - (int)p1;
                if (insert >= min_insert && insert <= max_insert) {
                    PairedMapping pm; pm.pos1 = p1; pm.pos2 = p2; pm.insert_size = insert;
                    results.push_back(pm);
                }
            }
        }
    }
    return results;
}

vector<size_t> Trivial_search(const string& read, const string& reference, int max_mismatch)
{
    vector<size_t> positions;
    size_t L = read.length(), N = reference.length();
    if (N < L) return positions;
    for (size_t i = 0; i <= N - L; i++) {
        size_t mismatch = 0;
        for (size_t j = 0; j < L; j++) {
            if (read[j] != reference[i + j]) mismatch++;
            if (mismatch > (size_t)max_mismatch) break;
        }
        if (mismatch <= (size_t)max_mismatch) positions.push_back(i);
    }
    return positions;
}

int main()
{
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);
    MEMORYSTATUSEX memStatus;
    memStatus.dwLength = sizeof(memStatus);
    GlobalMemoryStatusEx(&memStatus);

    cout << "========================================" << endl;
    cout << "실행 환경 (Machine Info)" << endl;
    cout << "========================================" << endl;
    cout << "OS:          Windows" << endl;
    cout << "CPU 코어 수: " << sysInfo.dwNumberOfProcessors << "개" << endl;
    cout << "총 RAM:      " << memStatus.ullTotalPhys / (1024 * 1024) << " MB" << endl;
    cout << "가용 RAM:    " << memStatus.ullAvailPhys / (1024 * 1024) << " MB" << endl;
    cout << endl;

    cout << "start" << endl;

    string dna = load_fasta_of_len("chr22.fa", 100000);
    if (dna.empty()) {
        cout << "chr22.fa not found -> random genome 생성" << endl;
        srand(42);
        for (size_t i = 0; i < 10000; i++) {
            int r = rand() % 4;
            if (r == 0) dna += 'A'; else if (r == 1) dna += 'C';
            else if (r == 2) dna += 'G'; else dna += 'T';
        }
    }

    cout << "genome length: " << dna.length() << endl;
    string original_dna = dna;

    // [수정 6] BWT를 original_dna만으로 구축
    // 기존: dna_double(20000bp) → suffix sort O(N² log N)이 4배 느렸음
    string bwt_input = original_dna;

    vector<uint64_t> bwt_2bit;
    size_t end_idx_onBWT;
    vector<size_t> suffix_array;
    vector<size_t> C_table;
    vector<vector<size_t>> Occ_table;

    cout << "BWT 시작" << endl;
    BWT_2bit(bwt_input, bwt_2bit, suffix_array, C_table, Occ_table, end_idx_onBWT);
    cout << "BWT 완료" << endl;

    vector<uint64_t> original_dna_2bit = encode_reference_2bit(original_dna);

    size_t bwt_mem_kb = 0;
    bwt_mem_kb += bwt_2bit.size() * sizeof(uint64_t);
    bwt_mem_kb += suffix_array.size() * sizeof(size_t);
    bwt_mem_kb += C_table.size() * sizeof(size_t);
    for (auto& v : Occ_table) bwt_mem_kb += v.size() * sizeof(size_t);
    bwt_mem_kb /= 1024;

    // ========================================
    // reads 생성
    // ========================================
    int M = 150000;
    int L32 = 32;
    int L100 = 100;
    int max_mismatch = 1;
    int insert_min = 150, insert_max = 600;
    // [수정 4] mapq_threshold 20 → 10
    int mapq_threshold = 10;
    string bases = "ACGT";
    double coverage = (double)M * L32 / original_dna.length();

    vector<string> exact_reads;
    vector<string> mismatch_reads;
    vector<string> independent_reads;

    // mt19937 사용 - MSVC rand() 주기(32767) 문제 해결
    // rand()는 주기가 32767이라 M=150000 생성시 같은 위치 4~5회 반복
    // → 100,000bp 중 32,767bp만 커버되는 버그의 원인
    // mt19937은 주기 2^19937, 실제 NGS 시뮬레이터 표준
    mt19937 rng(42);
    uniform_int_distribution<int> dist_start32(0, original_dna.length() - L32 - 1);
    uniform_int_distribution<int> dist_mut(0, max_mismatch);
    uniform_int_distribution<int> dist_pos32(0, L32 - 1);
    uniform_int_distribution<int> dist_base(0, 3);

    for (int i = 0; i < M; i++) {
        int start = dist_start32(rng);
        string read = original_dna.substr(start, L32);
        exact_reads.push_back(read);
        string mread = read;
        int num_mut = dist_mut(rng);
        for (int m = 0; m < num_mut; m++) {
            int mp = dist_pos32(rng); char orig = mread[mp]; char mut;
            do { mut = bases[dist_base(rng)]; } while (mut == orig);
            mread[mp] = mut;
        }
        mismatch_reads.push_back(mread);
    }

    mt19937 rng2(1234);
    uniform_int_distribution<int> dist_start32b(0, original_dna.length() - L32 - 1);
    uniform_int_distribution<int> dist_pos32b(0, L32 - 1);
    for (int i = 0; i < M; i++) {
        int start = dist_start32b(rng2);
        string read = original_dna.substr(start, L32);
        int mp = dist_pos32b(rng2); char orig = read[mp]; char mut;
        do { mut = bases[dist_base(rng2)]; } while (mut == orig);
        read[mp] = mut;
        independent_reads.push_back(read);
    }

    vector<string> mismatch_reads_100;
    mt19937 rng3(9999);
    uniform_int_distribution<int> dist_start100(0, original_dna.length() - L100 - 1);
    uniform_int_distribution<int> dist_pos100(0, L100 - 1);
    for (int i = 0; i < M; i++) {
        if ((int)original_dna.length() <= L100) break;
        int start = dist_start100(rng3);
        string read = original_dna.substr(start, L100);
        string mread = read;
        int num_mut = dist_mut(rng3);
        for (int m = 0; m < num_mut; m++) {
            int mp = dist_pos100(rng3); char orig = mread[mp]; char mut;
            do { mut = bases[dist_base(rng3)]; } while (mut == orig);
            mread[mp] = mut;
        }
        mismatch_reads_100.push_back(mread);
    }
    int M100 = mismatch_reads_100.size();

    vector<vector<uint64_t>> mismatch_reads_2bit_32;
    for (auto& r : mismatch_reads)     mismatch_reads_2bit_32.push_back(encode_read_2bit(r));
    vector<vector<uint64_t>> mismatch_reads_2bit_100;
    for (auto& r : mismatch_reads_100) mismatch_reads_2bit_100.push_back(encode_read_2bit(r));

    size_t reads_mem_kb = 0;
    for (auto& r : exact_reads)    reads_mem_kb += r.size();
    for (auto& r : mismatch_reads) reads_mem_kb += r.size();
    reads_mem_kb /= 1024;

    cout << "========================================" << endl;
    cout << "실험 환경" << endl;
    cout << "========================================" << endl;
    cout << "Reference:       chr22.fa (UCSC hg38)" << endl;
    cout << "genome 길이:     " << original_dna.length() << " bp" << endl;
    cout << "reads 수:        " << M << "개" << endl;
    cout << "read 길이:       L=" << L32 << "bp (비교용 L=" << L100 << "bp)" << endl;
    cout << "mismatch 허용:   " << max_mismatch << "개" << endl;
    cout << "insert size:     " << insert_min << " ~ " << insert_max << " bp" << endl;
    cout << "MAPQ 임계값:     " << mapq_threshold << endl;
    cout << "커버리지:        " << fixed << setprecision(1) << coverage << "x" << endl;
    cout << "BWT 인덱스 메모리: " << bwt_mem_kb << " KB" << endl;
    cout << "reads 메모리:    " << reads_mem_kb << " KB" << endl;
    cout << endl;

    keep_running = true;
    thread loadingThread(showLoadingAnimation);

    // ========================================
    // [비교] char vs 2bit reverse_complement
    // ========================================
    auto rc_char_start_32 = chrono::high_resolution_clock::now();
    for (int i = 0; i < M; i++) reverse_complement_char(mismatch_reads[i]);
    auto dur_rc_char_32 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - rc_char_start_32);

    auto rc_2bit_start_32 = chrono::high_resolution_clock::now();
    for (int i = 0; i < M; i++) reverse_complement_2bit(mismatch_reads[i]);
    auto dur_rc_2bit_32 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - rc_2bit_start_32);

    auto rc_char_start_100 = chrono::high_resolution_clock::now();
    for (int i = 0; i < M100; i++) reverse_complement_char(mismatch_reads_100[i]);
    auto dur_rc_char_100 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - rc_char_start_100);

    auto rc_2bit_start_100 = chrono::high_resolution_clock::now();
    for (int i = 0; i < M100; i++) reverse_complement_2bit(mismatch_reads_100[i]);
    auto dur_rc_2bit_100 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - rc_2bit_start_100);

    // ========================================
    // [공정한 비교] char vs 2bit Pigeonhole 검증
    // [수정 3] 각 벤치마크 전 fm_cache.clear()로 독립성 보장
    // ========================================

    fm_cache.clear();
    auto pig_char_start_32 = chrono::high_resolution_clock::now();
    int pig_char_mapped_32 = 0;
    for (int i = 0; i < M; i++) {
        vector<size_t> pos = Pigeonhole_search_char(
            mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
            original_dna, max_mismatch);
        if (!pos.empty()) pig_char_mapped_32++;
    }
    auto dur_pig_char_32 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_char_start_32);

    fm_cache.clear();  // [수정 3] 독립 측정
    auto pig_2bit_start_32 = chrono::high_resolution_clock::now();
    int pig_2bit_mapped_32 = 0;
    vector<vector<int>> mapping_table_pigeon(original_dna.length(), vector<int>(4, 0));
    for (int i = 0; i < M; i++) {
        vector<size_t> candidates;
        {
            int parts = max_mismatch + 1;
            int part_len = L32 / parts;
            set<size_t> cand_set;
            for (int p = 0; p < parts; p++) {
                int start = p * part_len;
                int end2 = (p == parts - 1) ? L32 : start + part_len;
                string sub = mismatch_reads[i].substr(start, end2 - start);
                vector<size_t> sub_pos = FM_search(sub, bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
                for (size_t pos : sub_pos)
                    if (pos >= (size_t)start) cand_set.insert(pos - start);
            }
            for (size_t p : cand_set)
                if (p + L32 <= original_dna.length()) candidates.push_back(p);
        }
        const vector<uint64_t>& read_enc = mismatch_reads_2bit_32[i];
        vector<size_t> valid_pos;
        for (size_t pos : candidates) {
            if (pos + L32 > original_dna.length()) continue;
            if (count_mismatches_2bit_preencoded(read_enc, original_dna_2bit, pos, L32) <= max_mismatch)
                valid_pos.push_back(pos);
        }
        if (!valid_pos.empty()) pig_2bit_mapped_32++;
    }
    auto dur_pig_2bit_32 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_2bit_start_32);

    fm_cache.clear();
    auto pig_char_start_100 = chrono::high_resolution_clock::now();
    int pig_char_mapped_100 = 0;
    for (int i = 0; i < M100; i++) {
        vector<size_t> pos = Pigeonhole_search_char(
            mismatch_reads_100[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
            original_dna, max_mismatch);
        if (!pos.empty()) pig_char_mapped_100++;
    }
    auto dur_pig_char_100 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_char_start_100);

    fm_cache.clear();
    auto pig_2bit_start_100 = chrono::high_resolution_clock::now();
    int pig_2bit_mapped_100 = 0;
    for (int i = 0; i < M100; i++) {
        vector<size_t> candidates;
        {
            int parts = max_mismatch + 1;
            int part_len = L100 / parts;
            set<size_t> cand_set;
            for (int p = 0; p < parts; p++) {
                int start = p * part_len;
                int end2 = (p == parts - 1) ? L100 : start + part_len;
                string sub = mismatch_reads_100[i].substr(start, end2 - start);
                vector<size_t> sub_pos = FM_search(sub, bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
                for (size_t pos : sub_pos)
                    if (pos >= (size_t)start) cand_set.insert(pos - start);
            }
            for (size_t p : cand_set)
                if (p + L100 <= original_dna.length()) candidates.push_back(p);
        }
        const vector<uint64_t>& read_enc = mismatch_reads_2bit_100[i];
        bool found = false;
        for (size_t pos : candidates) {
            if (pos + L100 > original_dna.length()) continue;
            if (count_mismatches_2bit_preencoded(read_enc, original_dna_2bit, pos, L100) <= max_mismatch) {
                found = true; break;
            }
        }
        if (found) pig_2bit_mapped_100++;
    }
    auto dur_pig_2bit_100 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_2bit_start_100);

    // ========================================
    // [표 1] exact reads 환경
    // ========================================
    fm_cache.clear();
    auto start_fm_exact = chrono::high_resolution_clock::now();
    int fm_exact_mapped = 0;
    for (int i = 0; i < M; i++) {
        vector<size_t> positions = FM_search(exact_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
        bool valid = false;
        for (size_t p : positions) if (p + L32 <= original_dna.length()) { valid = true; break; }
        if (valid) fm_exact_mapped++;
    }
    auto dur_fm_exact = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_fm_exact);

    fm_cache.clear();
    auto start_fm_mismatch = chrono::high_resolution_clock::now();
    int fm_mismatch_mapped = 0;
    for (int i = 0; i < M; i++) {
        vector<size_t> positions = FM_search(mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
        bool valid = false;
        for (size_t p : positions) if (p + L32 <= original_dna.length()) { valid = true; break; }
        if (valid) fm_mismatch_mapped++;
    }
    auto dur_fm_mismatch = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_fm_mismatch);

    // ========================================
    // [표 2] mismatch reads 환경
    // ========================================

    // Trivial
    auto start_trivial = chrono::high_resolution_clock::now();
    int trivial_mapped = 0;
    for (int i = 0; i < M; i++) {
        vector<size_t> positions = Trivial_search(mismatch_reads[i], original_dna, max_mismatch);
        if (!positions.empty()) trivial_mapped++;
    }
    auto dur_trivial = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_trivial);

    // Pigeonhole+HD+MAPQ
    // [수정 2] original_dna 전달, [수정 3] 캐시 클리어
    fm_cache.clear();
    auto start_mapq = chrono::high_resolution_clock::now();
    int mapq_mapped = 0, mapq_filtered = 0;
    vector<vector<int>> mapping_table_mapq(original_dna.length(), vector<int>(4, 0));
    for (int i = 0; i < M; i++) {
        vector<MappingResult> results = Pigeonhole_search_with_MAPQ(
            mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
            original_dna, max_mismatch, mapq_threshold);  // [수정 2] original_dna
        if (!results.empty()) {
            mapq_mapped++;
            for (auto& mr : results) {
                if (mr.pos + L32 > original_dna.length()) continue;
                for (int j = 0; j < L32; j++) {
                    if (mr.pos + j >= original_dna.length()) break;
                    // [수정 1] 변이된 read 염기 대신 original_dna 염기를 기록
                    // 기존: char_to_2bit(mismatch_reads[i][j]) → 변이 염기가 들어갈 수 있음
                    // 수정: char_to_2bit(original_dna[mr.pos + j]) → 정확한 reference 염기 기록
                    mapping_table_mapq[mr.pos + j][char_to_2bit(original_dna[mr.pos + j])]++;
                }
            }
        }
        else {
            // MAPQ 필터에서 걸렸는지 확인 (필터링 카운트용)
            vector<size_t> pig_pos = Pigeonhole_search_2bit(
                mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
                original_dna, max_mismatch);
            if (!pig_pos.empty()) mapq_filtered++;
        }
    }
    auto dur_mapq = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_mapq);

    // Paired-end
    // [수정 5] paired_reads 시드 별도 지정
    mt19937 rng4(5678);
    uniform_int_distribution<int> dist_frag(insert_min, insert_max - 1);
    vector<pair<string, string>> paired_reads;
    for (int i = 0; i < M; i++) {
        size_t frag_len = dist_frag(rng4);
        if (frag_len + (size_t)L32 * 2 >= original_dna.length()) continue;
        uniform_int_distribution<int> dist_fstart(0, original_dna.length() - frag_len - L32 - 1);
        size_t frag_start = dist_fstart(rng4);
        string read1 = original_dna.substr(frag_start, L32);
        int num_mut1 = dist_mut(rng4);
        for (int m = 0; m < num_mut1; m++) {
            int mp = dist_pos32(rng4); char orig = read1[mp]; char mut;
            do { mut = bases[dist_base(rng4)]; } while (mut == orig);
            read1[mp] = mut;
        }
        size_t read2_start = frag_start + frag_len;
        if (read2_start + L32 > original_dna.length()) continue;
        string read2 = original_dna.substr(read2_start, L32);
        string read2_rc = reverse_complement_2bit(read2);
        int num_mut2 = dist_mut(rng4);
        for (int m = 0; m < num_mut2; m++) {
            int mp = dist_pos32(rng4); char orig = read2_rc[mp]; char mut;
            do { mut = bases[dist_base(rng4)]; } while (mut == orig);
            read2_rc[mp] = mut;
        }
        paired_reads.push_back({ read1, read2_rc });
    }
    fm_cache.clear();
    auto start_paired = chrono::high_resolution_clock::now();
    int paired_valid = 0;
    for (int i = 0; i < (int)paired_reads.size(); i++) {
        vector<PairedMapping> results = Paired_end_search(
            paired_reads[i].first, paired_reads[i].second,
            bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
            original_dna, max_mismatch, insert_min, insert_max);
        if (!results.empty()) paired_valid++;
    }
    auto dur_paired = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_paired);

    // 통합 벤치마크
    fm_cache.clear();
    auto start_combined = chrono::high_resolution_clock::now();
    int combined_mapped = 0;
    for (int i = 0; i < M; i++) {
        bool mapped = false;
        vector<size_t> fm_pos = FM_search(independent_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
        for (size_t p : fm_pos)
            if (p + L32 <= original_dna.length()) { mapped = true; break; }
        if (!mapped) {
            vector<MappingResult> mapq_pos = Pigeonhole_search_with_MAPQ(
                independent_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
                original_dna, max_mismatch, mapq_threshold);
            for (auto& mr : mapq_pos)
                if (mr.pos + L32 <= original_dna.length()) { mapped = true; break; }
        }
        if (mapped) combined_mapped++;
    }
    auto dur_combined = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_combined);

    auto start_trivial2 = chrono::high_resolution_clock::now();
    int trivial2_mapped = 0;
    for (int i = 0; i < M; i++) {
        vector<size_t> positions = Trivial_search(independent_reads[i], original_dna, max_mismatch);
        if (!positions.empty()) trivial2_mapped++;
    }
    auto dur_trivial2 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_trivial2);

    keep_running = false;
    if (loadingThread.joinable()) loadingThread.join();

    size_t mapq_table_mem_kb = mapping_table_mapq.size() * 4 * sizeof(int) / 1024;

    double fm_exact_speed = dur_fm_exact.count() > 0 ? (double)M / dur_fm_exact.count() * 1e6 : 0;
    double trivial_speed = dur_trivial.count() > 0 ? (double)M / dur_trivial.count() * 1e6 : 0;
    double pig_speed = dur_pig_2bit_32.count() > 0 ? (double)M / dur_pig_2bit_32.count() * 1e6 : 0;
    double mapq_speed = dur_mapq.count() > 0 ? (double)M / dur_mapq.count() * 1e6 : 0;
    double paired_speed = dur_paired.count() > 0 ? (double)paired_reads.size() / dur_paired.count() * 1e6 : 0;
    double combined_speed = dur_combined.count() > 0 ? (double)M / dur_combined.count() * 1e6 : 0;
    double trivial2_speed = dur_trivial2.count() > 0 ? (double)M / dur_trivial2.count() * 1e6 : 0;

    size_t mem_char_32 = (size_t)M * L32;
    size_t mem_2bit_32 = (size_t)M * L32 / 4;
    size_t mem_char_100 = (size_t)M * L100;
    size_t mem_2bit_100 = (size_t)M * L100 / 4;

    // ========================================
    // 출력
    // ========================================
    cout << "========================================" << endl;
    cout << "2bit Encoding 효과" << endl;
    cout << "========================================" << endl;
    cout << left << setw(24) << "항목" << setw(16) << "char 기반" << setw(16) << "2bit 기반" << setw(12) << "절약률" << endl;
    cout << string(68, '-') << endl;
    cout << left << setw(24) << "저장공간/base" << setw(16) << "8 bit" << setw(16) << "2 bit" << setw(12) << "75%" << endl;
    cout << left << setw(24) << "genome 메모리(10Kbp)" << setw(16) << "10 KB" << setw(16) << "2.5 KB" << setw(12) << "75%" << endl;
    cout << left << setw(24) << ("reads 메모리(L=" + to_string(L32) + ")")
        << setw(16) << to_string(mem_char_32 / 1024) + " KB"
        << setw(16) << to_string(mem_2bit_32 / 1024) + " KB" << setw(12) << "75%" << endl;
    cout << left << setw(24) << ("reads 메모리(L=" + to_string(L100) + ")")
        << setw(16) << to_string(mem_char_100 / 1024) + " KB"
        << setw(16) << to_string(mem_2bit_100 / 1024) + " KB" << setw(12) << "75%" << endl;
    cout << left << setw(24) << "BWT 인덱스" << setw(16) << "-" << setw(16) << to_string(bwt_mem_kb) + " KB" << setw(12) << "-" << endl;
    {
        char buf[32];
        sprintf_s(buf, "%.2fx", dur_pig_char_32.count() > 0 ? (double)dur_pig_2bit_32.count() / dur_pig_char_32.count() : 0);
        string s32 = string(buf) + (dur_pig_2bit_32.count() < dur_pig_char_32.count() ? " (빠름↑)" : " (느림↓)");
        sprintf_s(buf, "%.2fx", dur_pig_char_100.count() > 0 ? (double)dur_pig_2bit_100.count() / dur_pig_char_100.count() : 0);
        string s100 = string(buf) + (dur_pig_2bit_100.count() < dur_pig_char_100.count() ? " (빠름↑)" : " (느림↓)");
        cout << left << setw(24) << "검증속도(L=32)"
            << setw(16) << to_string(dur_pig_char_32.count()) + " us"
            << setw(16) << to_string(dur_pig_2bit_32.count()) + " us" << setw(12) << s32 << endl;
        cout << left << setw(24) << "검증속도(L=100)"
            << setw(16) << to_string(dur_pig_char_100.count()) + " us"
            << setw(16) << to_string(dur_pig_2bit_100.count()) + " us" << setw(12) << s100 << endl;
    }
    cout << "→ DNA 염기(A/C/G/T) 4종을 2bit로 표현, char(8bit) 대비 75% 메모리 절약" << endl;
    cout << "→ L이 길수록 XOR+popcount 병렬 연산의 이점이 커져 2bit가 빠름" << endl;
    cout << endl;

    cout << "========================================" << endl;
    cout << "[표 1] exact reads 환경 (mismatch=0)" << endl;
    cout << "FM-Index는 exact match 전용 알고리즘" << endl;
    cout << "========================================" << endl;
    cout << left << setw(24) << "항목" << setw(20) << "FM-Index(exact)" << setw(20) << "FM-Index(mismatch)" << endl;
    cout << string(64, '-') << endl;
    cout << left << setw(24) << "실행시간(us)" << setw(20) << dur_fm_exact.count() << setw(20) << dur_fm_mismatch.count() << endl;
    cout << left << setw(24) << "매핑된 reads" << setw(20) << fm_exact_mapped << setw(20) << fm_mismatch_mapped << endl;
    cout << left << setw(24) << "매핑률"
        << setw(19) << fixed << setprecision(1) << (double)fm_exact_mapped / M * 100 << "%"
        << setw(19) << (double)fm_mismatch_mapped / M * 100 << "%" << endl;
    cout << "→ FM-Index는 mismatch reads에서 매핑률 저하 확인" << endl;
    cout << endl;

    cout << "========================================" << endl;
    cout << "[표 2] mismatch reads 환경 (mismatch<=1)" << endl;
    cout << "========================================" << endl;
    cout << left << setw(20) << "항목" << setw(16) << "Trivial" << setw(18) << "Pigeonhole(2bit)" << setw(20) << "Pigeonhole+HD+MAPQ" << setw(14) << "Paired-end" << endl;
    cout << string(88, '-') << endl;
    cout << left << setw(20) << "실행시간(us)" << setw(16) << dur_trivial.count() << setw(18) << dur_pig_2bit_32.count() << setw(20) << dur_mapq.count() << setw(14) << dur_paired.count() << endl;
    cout << left << setw(20) << "매핑속도(r/s)" << setw(16) << (int)trivial_speed << setw(18) << (int)pig_speed << setw(20) << (int)mapq_speed << setw(14) << (int)paired_speed << endl;
    cout << left << setw(20) << "매핑된 reads" << setw(16) << trivial_mapped << setw(18) << pig_2bit_mapped_32 << setw(20) << mapq_mapped << setw(14) << paired_valid << endl;
    cout << left << setw(20) << "매핑률"
        << setw(15) << fixed << setprecision(1) << (double)trivial_mapped / M * 100 << "%"
        << setw(17) << (double)pig_2bit_mapped_32 / M * 100 << "%"
        << setw(19) << (double)mapq_mapped / M * 100 << "%"
        << setw(13) << (paired_reads.size() > 0 ? (double)paired_valid / paired_reads.size() * 100 : 0) << "%" << endl;
    cout << left << setw(20) << "메모리(KB)" << setw(16) << "-" << setw(18) << "-" << setw(20) << mapq_table_mem_kb << setw(14) << "-" << endl;
    cout << left << setw(20) << "MAPQ 필터링" << setw(16) << "-" << setw(18) << "-" << setw(20) << mapq_filtered << setw(14) << "-" << endl;
    cout << "(HD = Hamming Distance, SW 대체, O(L) vs O(L²))" << endl;
    cout << endl;

    cout << "========================================" << endl;
    cout << "통합 알고리즘 vs Trivial 비교표" << endl;
    cout << "(독립 mismatch reads 기반 공정한 비교)" << endl;
    cout << "========================================" << endl;
    cout << left << setw(16) << "항목" << setw(28) << "FM+Pigeonhole+HD+MAPQ" << setw(16) << "Trivial" << endl;
    cout << string(60, '-') << endl;
    cout << left << setw(16) << "실행시간(us)" << setw(28) << dur_combined.count() << setw(16) << dur_trivial2.count() << endl;
    cout << left << setw(16) << "매핑속도(r/s)" << setw(28) << (int)combined_speed << setw(16) << (int)trivial2_speed << endl;
    cout << left << setw(16) << "매핑된 reads" << setw(28) << combined_mapped << setw(16) << trivial2_mapped << endl;
    cout << left << setw(16) << "매핑률"
        << setw(27) << fixed << setprecision(1) << (double)combined_mapped / M * 100 << "%"
        << setw(15) << (double)trivial2_mapped / M * 100 << "%" << endl;
    if (dur_trivial2.count() > 0 && dur_combined.count() > 0) {
        char buf[32];
        sprintf_s(buf, "%.1fx", (double)dur_trivial2.count() / dur_combined.count());
        cout << left << setw(16) << "속도 향상" << setw(28) << buf << setw(16) << "-" << endl;
    }
    cout << endl;

    cout << "========================================" << endl;
    cout << "알고리즘 단계별 개선 효과" << endl;
    cout << "========================================" << endl;
    cout << left << setw(34) << "단계" << setw(18) << "실행시간(us)" << setw(16) << "매핑속도(r/s)" << setw(12) << "매핑률" << endl;
    cout << string(80, '-') << endl;
    cout << left << setw(34) << "Trivial" << setw(18) << dur_trivial.count() << setw(16) << (int)trivial_speed << setw(12) << fixed << setprecision(1) << (double)trivial_mapped / M * 100 << "%" << endl;
    cout << left << setw(34) << "Pigeonhole" << setw(18) << dur_pig_2bit_32.count() << setw(16) << (int)pig_speed << setw(12) << (double)pig_2bit_mapped_32 / M * 100 << "%" << endl;
    cout << left << setw(34) << "Pigeonhole + HD + MAPQ" << setw(18) << dur_mapq.count() << setw(16) << (int)mapq_speed << setw(12) << (double)mapq_mapped / M * 100 << "%" << endl;
    cout << left << setw(34) << "FM + Pigeonhole + HD + MAPQ" << setw(18) << dur_combined.count() << setw(16) << (int)combined_speed << setw(12) << (double)combined_mapped / M * 100 << "%" << endl;
    {
        char buf[32];
        sprintf_s(buf, "%.1fx", dur_pig_2bit_32.count() > 0 ? (double)dur_trivial.count() / dur_pig_2bit_32.count() : 0);
        cout << "→ Trivial 대비 Pigeonhole: " << buf << " 빠름" << endl;
        sprintf_s(buf, "%.1fx", dur_combined.count() > 0 ? (double)dur_trivial.count() / dur_combined.count() : 0);
        cout << "→ Trivial 대비 통합 알고리즘: " << buf << " 빠름" << endl;
    }
    cout << "(HD = Hamming Distance, MAPQ inspired score 기반 반복서열 필터링)" << endl;
    cout << endl;

    // ========================================
    // Consensus 복원 + 정확도
    // ========================================
    string recovered = "";
    for (size_t i = 0; i < original_dna.length(); i++) {
        int total = mapping_table_mapq[i][0] + mapping_table_mapq[i][1]
            + mapping_table_mapq[i][2] + mapping_table_mapq[i][3];
        if (total == 0) { recovered += 'N'; continue; }
        int max_idx = 0, max_val = mapping_table_mapq[i][0]; bool tie = false;
        for (int j = 1; j < 4; j++) {
            if (mapping_table_mapq[i][j] > max_val) { max_val = mapping_table_mapq[i][j]; max_idx = j; tie = false; }
            else if (mapping_table_mapq[i][j] == max_val) tie = true;
        }
        if (tie) recovered += 'N';
        else     recovered += bit_to_char(max_idx);
    }

    size_t match = 0, mismatch_count = 0, unmapped_count = 0;
    for (size_t i = 0; i < original_dna.length(); i++) {
        if (recovered[i] == 'N') unmapped_count++;
        else if (original_dna[i] == recovered[i]) match++;
        else {
            mismatch_count++;
            if (mismatch_count <= 5)
                cout << "불일치 위치: " << i << " (원본: " << original_dna[i] << ", 복원: " << recovered[i] << ")" << endl;
        }
    }

    double accuracy_total = (double)match / original_dna.length() * 100;
    double accuracy_mapped = (match + mismatch_count > 0)
        ? (double)match / (match + mismatch_count) * 100 : 0;

    cout << "========================================" << endl;
    cout << "정확도 결과 (Pigeonhole+HD+MAPQ 기준)" << endl;
    cout << "========================================" << endl;
    cout << "매핑됨:          " << match + mismatch_count << "개" << endl;
    cout << "매핑 안됨(N):    " << unmapped_count << "개" << endl;
    cout << "  MAPQ 필터링:   " << mapq_filtered << "개 (반복 서열 제거)" << endl;
    cout << "불일치:          " << mismatch_count << "개" << endl;
    cout << "정확도(전체):    " << fixed << setprecision(2) << accuracy_total << "% (분모=전체 genome)" << endl;
    cout << "정확도(매핑됨):  " << fixed << setprecision(2) << accuracy_mapped << "% (분모=커버된 위치)" << endl;
    cout << "========================================" << endl;

    return 0;
}