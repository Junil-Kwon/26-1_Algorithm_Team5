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
#pragma comment(lib, "psapi.lib")

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
    for (uint8_t i = 0; i < 5; i++) { size_t count = C_table[i]; C_table[i] = total; total += count; }
}

vector<size_t> FM_search(const string& pattern, const vector<uint64_t>& bwt_packed,
    size_t end_idx_onBWT, const vector<size_t>& C_table,
    const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array)
{
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
        if (top >= bot) return positions;
    }
    for (size_t i = top; i < bot; i++) positions.push_back(suffix_array[i]);
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
    set<size_t> candidates;
    for (int p = 0; p < parts; p++) {
        int start = p * part_len;
        int end = (p == parts - 1) ? L : start + part_len;
        string sub = read.substr(start, end - start);
        vector<size_t> sub_positions = FM_search(sub, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array);
        for (size_t pos : sub_positions) {
            if (pos < (size_t)start) continue;
            candidates.insert(pos - start);
        }
    }
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
    set<size_t> candidates;
    for (int p = 0; p < parts; p++) {
        int start = p * part_len;
        int end = (p == parts - 1) ? L : start + part_len;
        string sub = read.substr(start, end - start);
        vector<size_t> sub_positions = FM_search(sub, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array);
        for (size_t pos : sub_positions) {
            if (pos < (size_t)start) continue;
            candidates.insert(pos - start);
        }
    }
    vector<uint64_t> read_encoded = encode_read_2bit(read);
    for (size_t pos : candidates) {
        if (pos + L > reference.length()) continue;
        if (count_mismatches_2bit(read_encoded, reference, pos, L) <= max_mismatch)
            positions.push_back(pos);
    }
    return positions;
}

struct MappingResult { size_t pos; int score; int mapq; };

vector<MappingResult> Pigeonhole_search_with_MAPQ(
    const string& read, const vector<uint64_t>& bwt_packed,
    size_t end_idx_onBWT, const vector<size_t>& C_table,
    const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array,
    const string& reference, int max_mismatch, int mapq_threshold = 20)
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
    const string& reference, const string& ref_double,
    int max_mismatch, int min_insert, int max_insert)
{
    vector<PairedMapping> results;
    size_t N = reference.length();
    vector<size_t> pos1_raw = Pigeonhole_search_2bit(read1, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array, ref_double, max_mismatch);
    vector<size_t> pos2_raw = Pigeonhole_search_2bit(read2_rc, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array, ref_double, max_mismatch);
    vector<size_t> pos1_list;
    for (size_t p : pos1_raw) if (p + read1.length() <= N) pos1_list.push_back(p);
    vector<size_t> pos2_list;
    string read2_original = reverse_complement_2bit(read2_rc);
    vector<uint64_t> read2_encoded = encode_read_2bit(read2_original);
    for (size_t p : pos2_raw) {
        if (p >= N && p + read2_rc.length() <= 2 * N) {
            size_t k = p - N;
            size_t L2 = read2_rc.length();
            if (k + L2 <= N) {
                size_t orig_pos = N - k - L2;
                if (count_mismatches_2bit(read2_encoded, reference, orig_pos, L2) <= max_mismatch)
                    pos2_list.push_back(orig_pos);
            }
        }
    }
    for (size_t p1 : pos1_list) {
        for (size_t p2 : pos2_list) {
            if (p2 >= p1) {
                int insert = (int)p2 + (int)read2_rc.length() - (int)p1;
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
    for (size_t i = 0; i <= N - L; i++) {
        size_t mismatch = 0;
        for (size_t j = 0; j < L; j++) {
            if (read[j] != reference[i + j]) mismatch++;
            if (mismatch > max_mismatch) break;
        }
        if (mismatch <= max_mismatch) positions.push_back(i);
    }
    return positions;
}

int main()
{
    // ========================================
    // 머신 정보 출력
    // ========================================
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

    string dna = load_fasta_of_len("C:/Workspace/PJ_alg/26-1_Algorithm_Team5/chr22.fa", 10000);
    if (dna.empty()) {
        cout << "chr22.fa not found -> random genome" << endl;
        srand((unsigned)time(0));
        for (size_t i = 0; i < 10000; i++) {
            int r = rand() % 4;
            if (r == 0) dna += 'A'; else if (r == 1) dna += 'C';
            else if (r == 2) dna += 'G'; else dna += 'T';
        }
    }

    cout << "genome length: " << dna.length() << endl;
    string original_dna = dna;
    string dna_rc = reverse_complement_2bit(original_dna);
    string dna_double = original_dna + dna_rc;

    vector<uint64_t> bwt_2bit;
    size_t end_idx_onBWT;
    vector<size_t> suffix_array;
    vector<size_t> C_table;
    vector<vector<size_t>> Occ_table;

    cout << "BWT 시작" << endl;
    BWT_2bit(dna_double, bwt_2bit, suffix_array, C_table, Occ_table, end_idx_onBWT);
    cout << "BWT 완료" << endl;

    // BWT 인덱스 메모리 직접 계산
    size_t bwt_mem_kb = 0;
    bwt_mem_kb += bwt_2bit.size() * sizeof(uint64_t);
    bwt_mem_kb += suffix_array.size() * sizeof(size_t);
    bwt_mem_kb += C_table.size() * sizeof(size_t);
    for (auto& v : Occ_table) bwt_mem_kb += v.size() * sizeof(size_t);
    bwt_mem_kb /= 1024;

    // ========================================
    // reads 생성
    // ========================================
    int M = 15000;
    int L32 = 32;
    int L100 = 100;
    int max_mismatch = 1;
    int insert_min = 100, insert_max = 500;
    int mapq_threshold = 20;
    string bases = "ACGT";
    double coverage = (double)M * L32 / original_dna.length();

    // exact reads (FM-Index용, mismatch 없음)
    vector<string> exact_reads;
    // mismatch reads (Trivial/Pigeonhole/MAPQ/Paired-end용)
    vector<string> mismatch_reads;
    // 통합 벤치마크용 독립 reads
    vector<string> independent_reads;

    srand(time(0));
    for (int i = 0; i < M; i++) {
        int start = rand() % (original_dna.length() - L32);
        string read = original_dna.substr(start, L32);
        exact_reads.push_back(read);

        string mread = read;
        int num_mut = rand() % (max_mismatch + 1);
        for (int m = 0; m < num_mut; m++) {
            int mp = rand() % L32; char orig = mread[mp]; char mut;
            do { mut = bases[rand() % 4]; } while (mut == orig);
            mread[mp] = mut;
        }
        mismatch_reads.push_back(mread);
    }
    srand(time(0) + 12345);
    for (int i = 0; i < M; i++) {
        int start = rand() % (original_dna.length() - L32);
        string read = original_dna.substr(start, L32);
        int mp = rand() % L32; char orig = read[mp]; char mut;
        do { mut = bases[rand() % 4]; } while (mut == orig);
        read[mp] = mut;
        independent_reads.push_back(read);
    }

    // L=100 reads
    vector<string> mismatch_reads_100;
    srand(time(0) + 99999);
    for (int i = 0; i < M; i++) {
        if ((int)original_dna.length() <= L100) break;
        int start = rand() % (original_dna.length() - L100);
        string read = original_dna.substr(start, L100);
        string mread = read;
        int num_mut = rand() % (max_mismatch + 1);
        for (int m = 0; m < num_mut; m++) {
            int mp = rand() % L100; char orig = mread[mp]; char mut;
            do { mut = bases[rand() % 4]; } while (mut == orig);
            mread[mp] = mut;
        }
        mismatch_reads_100.push_back(mread);
    }
    int M100 = mismatch_reads_100.size();

    // reads 메모리 직접 계산
    size_t reads_mem_kb = 0;
    for (auto& r : exact_reads)    reads_mem_kb += r.size();
    for (auto& r : mismatch_reads) reads_mem_kb += r.size();
    reads_mem_kb /= 1024;

    // ========================================
    // 실험 환경 출력
    // ========================================
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
    std::thread loadingThread(showLoadingAnimation);

    // ========================================
    // [비교] char vs 2bit: reverse_complement
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
    // [비교] char vs 2bit: Pigeonhole mismatch 검증
    // ========================================
    auto pig_char_start_32 = chrono::high_resolution_clock::now();
    int pig_char_mapped_32 = 0;
    for (int i = 0; i < M; i++) {
        vector<size_t> pos = Pigeonhole_search_char(mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array, dna_double, max_mismatch);
        for (size_t p : pos) if (p + L32 <= original_dna.length()) { pig_char_mapped_32++; break; }
    }
    auto dur_pig_char_32 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_char_start_32);

    auto pig_2bit_start_32 = chrono::high_resolution_clock::now();
    int pig_2bit_mapped_32 = 0;
    vector<vector<int>> mapping_table_pigeon(original_dna.length(), vector<int>(4, 0));
    for (int i = 0; i < M; i++) {
        vector<size_t> positions = Pigeonhole_search_2bit(mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array, dna_double, max_mismatch);
        vector<size_t> valid_pos;
        for (size_t p : positions) if (p + L32 <= original_dna.length()) valid_pos.push_back(p);
        if (!valid_pos.empty()) pig_2bit_mapped_32++;
        for (size_t pos : valid_pos)
            for (int j = 0; j < L32; j++) {
                if (pos + j >= original_dna.length()) break;
                mapping_table_pigeon[pos + j][char_to_2bit(mismatch_reads[i][j])]++;
            }
    }
    auto dur_pig_2bit_32 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_2bit_start_32);

    auto pig_char_start_100 = chrono::high_resolution_clock::now();
    int pig_char_mapped_100 = 0;
    for (int i = 0; i < M100; i++) {
        vector<size_t> pos = Pigeonhole_search_char(mismatch_reads_100[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array, dna_double, max_mismatch);
        for (size_t p : pos) if (p + L100 <= original_dna.length()) { pig_char_mapped_100++; break; }
    }
    auto dur_pig_char_100 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_char_start_100);

    auto pig_2bit_start_100 = chrono::high_resolution_clock::now();
    int pig_2bit_mapped_100 = 0;
    for (int i = 0; i < M100; i++) {
        vector<size_t> pos = Pigeonhole_search_2bit(mismatch_reads_100[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array, dna_double, max_mismatch);
        for (size_t p : pos) if (p + L100 <= original_dna.length()) { pig_2bit_mapped_100++; break; }
    }
    auto dur_pig_2bit_100 = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - pig_2bit_start_100);

    // ========================================
    // [표 1] exact reads 환경 (FM-Index 전용)
    // FM-Index는 exact match만 가능 → exact reads 사용
    // ========================================
    auto start_fm_exact = chrono::high_resolution_clock::now();
    int fm_exact_mapped = 0;
    vector<vector<size_t>> mapping_table_fm(original_dna.length(), vector<size_t>(4, 0));
    for (int i = 0; i < M; i++) {
        vector<size_t> positions = FM_search(exact_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
        vector<size_t> valid_pos;
        for (size_t p : positions) if (p + L32 <= original_dna.length()) valid_pos.push_back(p);
        if (!valid_pos.empty()) fm_exact_mapped++;
        for (size_t pos : valid_pos)
            for (int j = 0; j < L32; j++) {
                if (pos + j >= original_dna.length()) break;
                mapping_table_fm[pos + j][char_to_2bit(exact_reads[i][j])]++;
            }
    }
    auto dur_fm_exact = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_fm_exact);

    // FM-Index에 mismatch reads 넣었을 때 (한계 확인용)
    auto start_fm_mismatch = chrono::high_resolution_clock::now();
    int fm_mismatch_mapped = 0;
    for (int i = 0; i < M; i++) {
        vector<size_t> positions = FM_search(mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
        vector<size_t> valid_pos;
        for (size_t p : positions) if (p + L32 <= original_dna.length()) valid_pos.push_back(p);
        if (!valid_pos.empty()) fm_mismatch_mapped++;
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
    auto start_mapq = chrono::high_resolution_clock::now();
    int mapq_mapped = 0, mapq_filtered = 0;
    vector<vector<int>> mapping_table_mapq(original_dna.length(), vector<int>(4, 0));
    for (int i = 0; i < M; i++) {
        vector<MappingResult> results = Pigeonhole_search_with_MAPQ(
            mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
            dna_double, max_mismatch, mapq_threshold);
        vector<MappingResult> valid_results;
        for (auto& mr : results) if (mr.pos + L32 <= original_dna.length()) valid_results.push_back(mr);
        if (!valid_results.empty()) {
            mapq_mapped++;
            for (auto& mr : valid_results)
                for (int j = 0; j < L32; j++) {
                    if (mr.pos + j >= original_dna.length()) break;
                    mapping_table_mapq[mr.pos + j][char_to_2bit(mismatch_reads[i][j])]++;
                }
        }
        else {
            vector<size_t> pig_pos = Pigeonhole_search_2bit(mismatch_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array, dna_double, max_mismatch);
            bool found = false;
            for (size_t p : pig_pos) if (p + L32 <= original_dna.length()) { found = true; break; }
            if (found) mapq_filtered++;
        }
    }
    auto dur_mapq = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_mapq);

    // Paired-end
    auto start_paired = chrono::high_resolution_clock::now();
    vector<pair<string, string>> paired_reads;
    srand(time(0));
    for (int i = 0; i < M; i++) {
        size_t frag_len = insert_min + rand() % (insert_max - insert_min);
        if (frag_len + L32 >= original_dna.length()) continue;
        size_t frag_start = rand() % (original_dna.length() - frag_len - L32);
        string read1 = original_dna.substr(frag_start, L32);
        int num_mut1 = rand() % (max_mismatch + 1);
        for (int m = 0; m < num_mut1; m++) {
            int mp = rand() % L32; char orig = read1[mp]; char mut;
            do { mut = bases[rand() % 4]; } while (mut == orig);
            read1[mp] = mut;
        }
        string read2 = original_dna.substr(frag_start + frag_len, L32);
        string read2_rc = reverse_complement_2bit(read2);
        int num_mut2 = rand() % (max_mismatch + 1);
        for (int m = 0; m < num_mut2; m++) {
            int mp = rand() % L32; char orig = read2_rc[mp]; char mut;
            do { mut = bases[rand() % 4]; } while (mut == orig);
            read2_rc[mp] = mut;
        }
        paired_reads.push_back({ read1, read2_rc });
    }
    int paired_valid = 0;
    for (int i = 0; i < (int)paired_reads.size(); i++) {
        vector<PairedMapping> results = Paired_end_search(
            paired_reads[i].first, paired_reads[i].second,
            bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
            original_dna, dna_double, max_mismatch, insert_min, insert_max);
        if (!results.empty()) paired_valid++;
    }
    auto dur_paired = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_paired);

    // ========================================
    // 통합 벤치마크 (FM + Pigeonhole+HD+MAPQ vs Trivial)
    // 독립 reads 기반 공정한 비교
    // ========================================
    auto start_combined = chrono::high_resolution_clock::now();
    int combined_mapped = 0;
    for (int i = 0; i < M; i++) {
        bool mapped = false;
        set<size_t> combined_positions;
        vector<size_t> fm_pos = FM_search(independent_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
        for (size_t p : fm_pos)
            if (p + L32 <= original_dna.length()) { combined_positions.insert(p); mapped = true; }
        vector<MappingResult> mapq_pos = Pigeonhole_search_with_MAPQ(
            independent_reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array,
            dna_double, max_mismatch, mapq_threshold);
        for (auto& mr : mapq_pos)
            if (mr.pos + L32 <= original_dna.length()) { combined_positions.insert(mr.pos); mapped = true; }
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

    // mapping_table 메모리 계산
    size_t mapq_table_mem_kb = mapping_table_mapq.size() * 4 * sizeof(int) / 1024;

    double fm_exact_speed = dur_fm_exact.count() > 0 ? (double)M / dur_fm_exact.count() * 1e6 : 0;
    double trivial_speed = dur_trivial.count() > 0 ? (double)M / dur_trivial.count() * 1e6 : 0;
    double pig_speed = dur_pig_2bit_32.count() > 0 ? (double)M / dur_pig_2bit_32.count() * 1e6 : 0;
    double mapq_speed = dur_mapq.count() > 0 ? (double)M / dur_mapq.count() * 1e6 : 0;
    double paired_speed = dur_paired.count() > 0 ? (double)M / dur_paired.count() * 1e6 : 0;
    double combined_speed = dur_combined.count() > 0 ? (double)M / dur_combined.count() * 1e6 : 0;
    double trivial2_speed = dur_trivial2.count() > 0 ? (double)M / dur_trivial2.count() * 1e6 : 0;

    // ========================================
    // 출력 1: char vs 2bit 비교
    // ========================================
    // char vs 2bit 메모리 계산
    // char: 1 문자 = 8bit = 1 byte
    // 2bit: 1 문자 = 2bit = 0.25 byte
    size_t mem_char_32 = (size_t)M * L32;          // byte
    size_t mem_2bit_32 = (size_t)M * L32 / 4;      // byte
    size_t mem_char_100 = (size_t)M * L100;
    size_t mem_2bit_100 = (size_t)M * L100 / 4;

    cout << "========================================" << endl;
    cout << "char vs 2bit: 메모리 & 속도 비교" << endl;
    cout << "========================================" << endl;
    cout << left
        << setw(10) << "L"
        << setw(16) << "char 속도(us)"
        << setw(16) << "2bit 속도(us)"
        << setw(18) << "속도 비율"
        << setw(16) << "char 메모리"
        << setw(16) << "2bit 메모리"
        << setw(12) << "절약률" << endl;
    cout << string(104, '-') << endl;
    {
        char buf[32];
        sprintf_s(buf, "%.2fx", dur_rc_char_32.count() > 0 ? (double)dur_rc_2bit_32.count() / dur_rc_char_32.count() : 0);
        cout << left << setw(10) << "L=32"
            << setw(16) << dur_rc_char_32.count()
            << setw(16) << dur_rc_2bit_32.count()
            << setw(18) << buf
            << setw(16) << to_string(mem_char_32 / 1024) + " KB"
            << setw(16) << to_string(mem_2bit_32 / 1024) + " KB"
            << setw(12) << "75%" << endl;
        sprintf_s(buf, "%.2fx", dur_rc_char_100.count() > 0 ? (double)dur_rc_2bit_100.count() / dur_rc_char_100.count() : 0);
        cout << left << setw(10) << "L=100"
            << setw(16) << dur_rc_char_100.count()
            << setw(16) << dur_rc_2bit_100.count()
            << setw(18) << buf
            << setw(16) << to_string(mem_char_100 / 1024) + " KB"
            << setw(16) << to_string(mem_2bit_100 / 1024) + " KB"
            << setw(12) << "75%" << endl;
    }
    cout << "※ 2bit: 속도↓ 메모리↑ (매 연산 인코딩 오버헤드 발생)" << endl;
    cout << "※ reference 자체를 2bit 저장 시 인코딩 비용 제거 가능" << endl;
    cout << endl;

    cout << "========================================" << endl;
    cout << "char vs 2bit: Pigeonhole mismatch 검증 비교" << endl;
    cout << "========================================" << endl;
    cout << left
        << setw(10) << "L"
        << setw(16) << "char(us)"
        << setw(16) << "2bit(us)"
        << setw(18) << "속도 비율"
        << setw(16) << "char 메모리"
        << setw(16) << "2bit 메모리"
        << setw(10) << "매핑률" << endl;
    cout << string(102, '-') << endl;
    {
        char buf[32];
        sprintf_s(buf, "%.2fx", dur_pig_char_32.count() > 0 ? (double)dur_pig_2bit_32.count() / dur_pig_char_32.count() : 0);
        cout << left << setw(10) << "L=32"
            << setw(16) << dur_pig_char_32.count()
            << setw(16) << dur_pig_2bit_32.count()
            << setw(18) << buf
            << setw(16) << to_string(mem_char_32 / 1024) + " KB"
            << setw(16) << to_string(mem_2bit_32 / 1024) + " KB"
            << setw(10) << fixed << setprecision(1) << (double)pig_2bit_mapped_32 / M * 100 << "%" << endl;
        sprintf_s(buf, "%.2fx", dur_pig_char_100.count() > 0 ? (double)dur_pig_2bit_100.count() / dur_pig_char_100.count() : 0);
        cout << left << setw(10) << "L=100"
            << setw(16) << dur_pig_char_100.count()
            << setw(16) << dur_pig_2bit_100.count()
            << setw(18) << buf
            << setw(16) << to_string(mem_char_100 / 1024) + " KB"
            << setw(16) << to_string(mem_2bit_100 / 1024) + " KB"
            << setw(10) << (double)pig_2bit_mapped_100 / M100 * 100 << "%" << endl;
    }
    cout << endl;

    // ========================================
    // 출력 2: [표 1] exact reads 환경
    // FM-Index는 exact match 전용임을 명확히 보여줌
    // ========================================
    cout << "========================================" << endl;
    cout << "[표 1] exact reads 환경 (mismatch=0)" << endl;
    cout << "FM-Index는 exact match 전용 알고리즘" << endl;
    cout << "========================================" << endl;
    cout << left
        << setw(24) << "항목"
        << setw(20) << "FM-Index(exact)"
        << setw(20) << "FM-Index(mismatch)" << endl;
    cout << string(64, '-') << endl;
    cout << left
        << setw(24) << "실행시간(us)"
        << setw(20) << dur_fm_exact.count()
        << setw(20) << dur_fm_mismatch.count() << endl;
    cout << left
        << setw(24) << "매핑속도(r/s)"
        << setw(20) << (int)fm_exact_speed
        << setw(20) << (int)(dur_fm_mismatch.count() > 0 ? (double)M / dur_fm_mismatch.count() * 1e6 : 0) << endl;
    cout << left
        << setw(24) << "매핑된 reads"
        << setw(20) << fm_exact_mapped
        << setw(20) << fm_mismatch_mapped << endl;
    cout << left
        << setw(24) << "매핑률"
        << setw(19) << fixed << setprecision(1) << (double)fm_exact_mapped / M * 100 << "%"
        << setw(19) << (double)fm_mismatch_mapped / M * 100 << "%" << endl;
    cout << "→ FM-Index는 mismatch reads에서 매핑률 저하 확인" << endl;
    cout << endl;

    // ========================================
    // 출력 3: [표 2] mismatch reads 환경
    // ========================================
    cout << "========================================" << endl;
    cout << "[표 2] mismatch reads 환경 (mismatch<=1)" << endl;
    cout << "========================================" << endl;
    cout << left
        << setw(20) << "항목"
        << setw(16) << "Trivial"
        << setw(18) << "Pigeonhole(2bit)"
        << setw(20) << "Pigeonhole+HD+MAPQ"
        << setw(14) << "Paired-end" << endl;
    cout << string(88, '-') << endl;
    cout << left
        << setw(20) << "실행시간(us)"
        << setw(16) << dur_trivial.count()
        << setw(18) << dur_pig_2bit_32.count()
        << setw(20) << dur_mapq.count()
        << setw(14) << dur_paired.count() << endl;
    cout << left
        << setw(20) << "매핑속도(r/s)"
        << setw(16) << (int)trivial_speed
        << setw(18) << (int)pig_speed
        << setw(20) << (int)mapq_speed
        << setw(14) << (int)paired_speed << endl;
    cout << left
        << setw(20) << "매핑된 reads"
        << setw(16) << trivial_mapped
        << setw(18) << pig_2bit_mapped_32
        << setw(20) << mapq_mapped
        << setw(14) << paired_valid << endl;
    cout << left
        << setw(20) << "매핑률"
        << setw(15) << fixed << setprecision(1) << (double)trivial_mapped / M * 100 << "%"
        << setw(17) << (double)pig_2bit_mapped_32 / M * 100 << "%"
        << setw(19) << (double)mapq_mapped / M * 100 << "%"
        << setw(13) << (double)paired_valid / (int)paired_reads.size() * 100 << "%" << endl;
    cout << left
        << setw(20) << "메모리(KB)"
        << setw(16) << "-"
        << setw(18) << "-"
        << setw(20) << mapq_table_mem_kb
        << setw(14) << "-" << endl;
    cout << left
        << setw(20) << "MAPQ 필터링"
        << setw(16) << "-"
        << setw(18) << "-"
        << setw(20) << mapq_filtered
        << setw(14) << "-" << endl;
    cout << "(HD = Hamming Distance, SW 대체, O(L) vs O(L²))" << endl;
    cout << endl;

    // ========================================
    // 출력 4: 통합 vs Trivial
    // ========================================
    cout << "========================================" << endl;
    cout << "통합 알고리즘 vs Trivial 비교표" << endl;
    cout << "(독립 mismatch reads 기반 공정한 비교)" << endl;
    cout << "========================================" << endl;
    cout << left
        << setw(16) << "항목"
        << setw(28) << "FM+Pigeonhole+HD+MAPQ"
        << setw(16) << "Trivial" << endl;
    cout << string(60, '-') << endl;
    cout << left
        << setw(16) << "실행시간(us)"
        << setw(28) << dur_combined.count()
        << setw(16) << dur_trivial2.count() << endl;
    cout << left
        << setw(16) << "매핑속도(r/s)"
        << setw(28) << (int)combined_speed
        << setw(16) << (int)trivial2_speed << endl;
    cout << left
        << setw(16) << "매핑된 reads"
        << setw(28) << combined_mapped
        << setw(16) << trivial2_mapped << endl;
    cout << left
        << setw(16) << "매핑률"
        << setw(27) << fixed << setprecision(1) << (double)combined_mapped / M * 100 << "%"
        << setw(15) << (double)trivial2_mapped / M * 100 << "%" << endl;
    if (dur_trivial2.count() > 0) {
        char buf[32];
        sprintf_s(buf, "%.1fx", (double)dur_trivial2.count() / dur_combined.count());
        cout << left << setw(16) << "속도 향상" << setw(28) << buf << setw(16) << "-" << endl;
    }
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
            if (mismatch_count <= 10)
                cout << "불일치 위치: " << i
                << " (원본: " << original_dna[i]
                << ", 복원: " << recovered[i] << ")" << endl;
        }
    }

    double accuracy = (double)match / original_dna.length() * 100;

    cout << "========================================" << endl;
    cout << "정확도 결과 (Pigeonhole+HD+MAPQ 기준)" << endl;
    cout << "========================================" << endl;
    cout << "매핑됨:          " << match + mismatch_count << "개" << endl;
    cout << "매핑 안됨(N):    " << unmapped_count << "개" << endl;
    cout << "  MAPQ 필터링:   " << mapq_filtered << "개 (반복 서열 제거)" << endl;
    cout << "불일치:          " << mismatch_count << "개" << endl;
    cout << "정확도:          " << fixed << setprecision(2) << accuracy << "%" << endl;
    cout << endl;
    cout << "[정확도 시각화]" << endl;
    cout << "0%    25%    50%    75%   100%" << endl;
    cout << "|-----|------|------|-----|" << endl;
    int bar_len = 25;
    int filled = (int)(accuracy / 100.0 * bar_len);
    cout << "|";
    for (int i = 0; i < bar_len; i++) cout << (i < filled ? "#" : " ");
    cout << "| " << fixed << setprecision(1) << accuracy << "%" << endl;
    cout << "========================================" << endl;

    return 0;
}