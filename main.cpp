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

using namespace std;
// 64bp 단위로 Occ 테이블 체크포인트를 설정 (64비트 정수 2개 단위)
const size_t CHECKPOINT_INTERVAL = 64;

string load_fasta_of_len(const string& filename, size_t max_len)
{
    ifstream file(filename); // FASTA 파일 열기

    if (!file.is_open())
    {
        cerr << "파일을 열 수 없습니다!" << endl;
        return "";
    }

    string line, genome_Seq;
    // 메모리 조각화를 방지하기 위해 미리 메모리 공간 예약
    genome_Seq.reserve(51000000); // chromosome 22 (약 5.08천만 bp)
    //genome_Seq.reserve(3200000000); // GRCh38 전체 유전체 (약 3.1억 bp)

    size_t str_len = 0;

    while (getline(file, line)) // 파일에서 한 줄씩 읽기
    {
        if (str_len >= max_len) break;

        if (line.empty()) continue;
        if (line[0] == '>')
        {
            // 헤더 라인(파일 설명)
            cout << "읽는 중: " << line << endl;
            continue;
        }

        for (char c : line)
        {
            c = toupper(c); // 소문자 대응
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
            {
                genome_Seq += c;
                str_len++;
            }
            // 'N'은 무시하고 저장하지 않음
        }
    }

    file.close();
    cout << endl << "성공적으로 파일을 읽어왔습니다." << endl << endl;

    return genome_Seq;
}

unsigned char char_to_2bit(char base)
{
    switch (base)
    {
    case 'A': return 0; // 00
    case 'C': return 1; // 01
    case 'G': return 2; // 10
    case 'T': return 3; // 11
        //case '$': return 0; // 끝 표시 문자는 00으로 처리 (특수 인덱스로 구분)
    default:  return 0; // 예외 처리
    }
}


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

//비트 연산 기반으로 특정 구간 내의 염기 개수를 카운팅하는 함수
size_t count_bits_in_chunk(uint64_t chunk, size_t active_bases, uint8_t target_2bit) {
    if (active_bases == 0) return 0;

    // 조사할 범위(active_bases)만큼만 상위 비트에서 마스킹하여 추출합니다.
    // 32개 미만일 때 오른쪽 쓰레기 비트들을 완전히 청소합니다.
    if (active_bases < 32) {
        size_t shift_amount = 64 - (active_bases * 2);
        chunk = (chunk >> shift_amount) << shift_amount;
    }

    uint64_t match_bits = 0;

    if (target_2bit == 0) { // 'A' (00)
        uint64_t inverted = ~chunk;
        match_bits = (inverted >> 1) & inverted & 0x5555555555555555ULL;
    }
    else if (target_2bit == 1) { // 'C' (01)
        match_bits = (~chunk >> 1) & chunk & 0x5555555555555555ULL;
    }
    else if (target_2bit == 2) { // 'G' (10)
        match_bits = (chunk >> 1) & ~chunk & 0x5555555555555555ULL;
    }
    else if (target_2bit == 3) { // 'T' (11)
        match_bits = (chunk >> 1) & chunk & 0x5555555555555555ULL;
    }

    // 최종 매칭된 비트들 중에서도 active_bases 범위 안의 것만 카운트하도록 유효 마스크 적용
    if (active_bases < 32) {
        uint64_t mask = ~((1ULL << (64 - active_bases * 2)) - 1);
        match_bits &= mask;
    }

    return __popcnt64(match_bits);
}

// [핵심 변경] 비트 배열과 체크포인트를 활용한 정밀 고속 Rank 계산 함수
size_t get_rank(size_t idx, uint8_t char_idx, const vector<uint64_t>& bwt_packed, const vector<vector<size_t>>& Occ_table, size_t end_idx_onBWT) {
    if (idx == 0) return 0;

    // $ 문자(char_idx == 0)의 고속 Rank 추적 보정
    if (char_idx == 0) {
        return (idx > end_idx_onBWT) ? 1 : 0;
    }

    // 나눗셈을 통해 포함되는 기저 블록과 남은 개수를 깔끔하게 분리
    size_t checkpoint_block = idx / CHECKPOINT_INTERVAL;
    size_t remaining = idx % CHECKPOINT_INTERVAL;

    // 기저 체크포인트 누적값 획득
    size_t rank_count = Occ_table[char_idx][checkpoint_block];

    // 만약 딱 64의 배수로 떨어져서 남은 구간이 없다면 그대로 반환
    if (remaining == 0) return rank_count;

    // 남은 구간이 있다면 해당 블록의 packed_idx부터 비트를 추적
    size_t packed_idx = checkpoint_block * 2;
    uint8_t target_2bit = char_idx - 1; // A=0, C=1, G=2, T=3

    if (remaining <= 32) {
        rank_count += count_bits_in_chunk(bwt_packed[packed_idx], remaining, target_2bit);
    }
    else {
        rank_count += count_bits_in_chunk(bwt_packed[packed_idx], 32, target_2bit);
        rank_count += count_bits_in_chunk(bwt_packed[packed_idx + 1], remaining - 32, target_2bit);
    }

    if (target_2bit == 0 && idx > end_idx_onBWT) {
        // 단, 기저 체크포인트(Occ_table)에는 이미 '$'를 분리해서 올바른 'A' 개수만 저장했으므로,
        // 현재 체크포인트 블록 '이후'에 '$'가 존재하면서 여전히 idx보다는 앞에 있을 때만 뺀다.
        size_t checkpoint_limit = checkpoint_block * CHECKPOINT_INTERVAL;
        if (end_idx_onBWT >= checkpoint_limit) {
            rank_count--;
        }
    }

    return rank_count;
}

// 2-bit BWT 수행 및 가소화용 Occ 테이블 구축
void BWT_2bit(string& text, vector<uint64_t>& bwt_packed, vector<size_t>& suffix_array, vector<size_t>& C_table, vector<vector<size_t>>& Occ_table, size_t& end_idx_onBWT) {
    text += '$';
    size_t length = text.length();

    suffix_array.clear();
    suffix_array.reserve(length);
    for (size_t i = 0; i < length; i++) suffix_array.push_back(i);

    sort(suffix_array.begin(), suffix_array.end(), [&](size_t a, size_t b) {
        for (size_t i = 0; i < length; i++) {
            char char_a = text[(a + i) % length];
            char char_b = text[(b + i) % length];
            if (char_a != char_b) return char_a < char_b;
        }
        return false;
        });

    C_table.assign(5, 0);

    // Occ_table의 가로 크기를 length가 아닌 (length / 64 + 1)로 대폭 압축 (메모리 절약)
    size_t num_checkpoints = length / CHECKPOINT_INTERVAL + 1;
    Occ_table.assign(5, vector<size_t>(num_checkpoints, 0));
    vector<size_t> current_occ(5, 0);

    bwt_packed.clear();
    bwt_packed.reserve(length / 32 + 1); // uint64_t 하나당 32개 염기 저장

    uint64_t buffer = 0;
    for (size_t i = 0; i < length; i++) {
        size_t last_idx = (suffix_array[i] + length - 1) % length;
        char suffix = text[last_idx];

        // $ = 0, A = 1, C = 2, G = 3, T = 4
        uint8_t char_idx = (suffix == '$') ? 0 : char_to_2bit(suffix) + 1;

        // 64bp 주기(체크포인트 블록의 시작점)마다 당시의 누적 빈도수를 마킹하여 저장
        if (i % CHECKPOINT_INTERVAL == 0) {
            size_t block = i / CHECKPOINT_INTERVAL;
            for (uint8_t b = 0; b < 5; b++) Occ_table[b][block] = current_occ[b];
        }

        current_occ[char_idx]++;
        C_table[char_idx]++;

        if (suffix == '$') end_idx_onBWT = i;

        // $ 문자는 비트 배열에서 00처리 되더라도, char_idx=0으로 Occ테이블에서 엄격히 격리되므로
        // 실제 염기(A,C,G,T)를 패킹할 때만 target_2bit를 올바르게 주입합니다.
        uint64_t b_val = (suffix == '$') ? 0 : char_to_2bit(suffix);
        size_t bit_offset = (i % 32) * 2;
        buffer |= (b_val << (62 - bit_offset));

        if (i % 32 == 32 - 1 || i == length - 1) {
            bwt_packed.push_back(buffer);
            buffer = 0;
        }
    }

    size_t total = 0;
    for (uint8_t i = 0; i < 5; i++) {
        size_t count = C_table[i];
        C_table[i] = total;
        total += count;
    }
}

// 비트와이즈 가속화 랭크 함수를 적용한 FM_search
vector<size_t> FM_search(const string& pattern, const vector<uint64_t>& bwt_packed, size_t end_idx_onBWT, const vector<size_t>& C_table, const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array) {
    vector<size_t> positions;
    size_t length = suffix_array.size();

    size_t top = 0;
    size_t bot = length;

    for (int i = pattern.length() - 1; i >= 0; i--) {
        char c = pattern[i];
        if (c == '$') return positions; // 패턴에 $가 들어오는 예외 차단
        uint8_t char_idx = char_to_2bit(c) + 1;

        // 고속 비트와이즈 랭크 계산
        size_t occ_top = get_rank(top, char_idx, bwt_packed, Occ_table, end_idx_onBWT);
        size_t occ_bot = get_rank(bot, char_idx, bwt_packed, Occ_table, end_idx_onBWT);

        top = C_table[char_idx] + occ_top;
        bot = C_table[char_idx] + occ_bot;

        if (top >= bot) return positions;
    }

    for (size_t i = top; i < bot; i++) {
        positions.push_back(suffix_array[i]);
    }
    return positions;
}

/*
// 2bit BWT에서 i번째 문자를 꺼내는 함수
char get_bwt_char(const vector<uint8_t>& bwt_2bit, size_t i, size_t end_idx)
{
    if (i == end_idx) return '$';

    size_t byte_idx = i / 4;
    size_t bit_offset = i % 4;

    uint8_t bits = (bwt_2bit[byte_idx] >> (6 - bit_offset * 2)) & 0x03;

    return bit_to_char(bits);
}
*/

/*
void BWT(string& text, string& bwt, vector<size_t>& suffix_array, vector<size_t>& C_table, vector<vector<size_t>>& Occ_table)
{
    text += '$'; // 끝 표시 문자 추가
    size_t length = text.length(); // 텍스트 길이 계산

    suffix_array.clear();
    suffix_array.reserve(length); // SA 메모리 예약

    for (size_t i = 0; i < length; i++) suffix_array.push_back(i); // SA 초기화 (0부터 length-1까지의 인덱스)

    // SA 배열 정렬 (접미사들을 기준으로 하여, 사전순으로 비교하여 정렬)
    //!!!NOTICE!!!: SA-IS 알고리즘을 사용하여, O(n^2logn) -> O(n)로 개선 가능, 구현 많이 복잡
    sort(suffix_array.begin(), suffix_array.end(), [&](size_t a, size_t b) {
        for (size_t i = 0; i < length; i++)
        {
            char char_a = text[(a + i) % length];
            char char_b = text[(b + i) % length];
            if (char_a != char_b) {
                return char_a < char_b;
            }
        }
        return false;
        });

    // C Table 준비
    C_table.assign(5, 0);

    // Occ Table 준비
    Occ_table.assign(5, vector<size_t>(length, 0));
    vector<size_t> current_occ(5, 0);

    bwt = ""; // BWT 문자열 초기화
    bwt.reserve(length); // 메모리 예약

    // BWT 문자열 생성
    for (size_t i = 0; i < length; i++)
    {
        size_t last_idx = (suffix_array[i] + length - 1) % length; // LF 매핑을 이용하여 마지막 인덱스 계산
        char suffix = text[last_idx]; // 접미사 추출

        uint8_t char_idx; // Occ 및 C 테이블에서 사용할 인덱스

        if (suffix == '$')
        {
            char_idx = 0; // $는 0번 인덱스
        }
        else
        {
            // A=1, C=2, G=3, T=4 (char_to_2bit 결과에 1을 더함)
            char_idx = char_to_2bit(suffix) + 1;
        }

        current_occ[char_idx]++;
        C_table[char_idx]++;

        // 현재 위치 i에 대한 Occ 값 기록
        for (int b = 0; b < 5; b++) Occ_table[b][i] = current_occ[b];

        bwt += text[last_idx]; // BWT 문자열에 접미사에서 추출한 문자를 추가
    }

    // C-Table 완성 (개수 -> 시작 인덱스로 변환)
    size_t total = 0;
    for (uint8_t i = 0; i < 5; i++)
    {
        size_t count = C_table[i];
        C_table[i] = total;
        total += count;
    }
}
*/

/*
void BWT_2bit(string& text, vector<uint8_t>& bwt_2bit, vector<size_t>& suffix_array, vector<size_t>& C_table, vector<vector<size_t>>& Occ_table, size_t& end_idx_onBWT)
{
    text += '$'; // 끝 표시 문자 추가
    size_t length = text.length(); // 텍스트 길이 계산

    suffix_array.clear();
    suffix_array.reserve(length); // SA 메모리 예약

    for (size_t i = 0; i < length; i++) suffix_array.push_back(i); // SA 초기화 (0부터 length-1까지의 인덱스)

    // SA 배열 정렬 (접미사들을 기준으로 하여, 사전순으로 비교하여 정렬)
    //!!!NOTICE!!!: SA-IS 알고리즘을 사용하여, O(n^2logn) -> O(n)로 개선 가능, 구현 많이 복잡
    sort(suffix_array.begin(), suffix_array.end(), [&](size_t a, size_t b) {
        for (size_t i = 0; i < length; i++)
        {
            char char_a = text[(a + i) % length];
            char char_b = text[(b + i) % length];
            if (char_a != char_b)
            {
                return char_a < char_b;
            }
        }
        return false;
        });

    // C Table 준비
    C_table.assign(5, 0);

    // Occ Table 준비
    // Occ Table은 체크포인트를 지정하여 저장하면 메모리를 절약할 수 있음.
    Occ_table.assign(5, vector<size_t>(length, 0));
    vector<size_t> current_occ(5, 0);

    bwt_2bit.clear();
    bwt_2bit.reserve(length / 4 + 1); // 4개당 1바이트

    unsigned char buffer = 0;
    bool special_char_found = false; // 끝 표시 문자를 찾았는지 여부
    for (size_t i = 0; i < length; i++)
    {
        size_t last_idx = (suffix_array[i] + length - 1) % length; // LF 매핑을 이용하여 마지막 인덱스 계산
        char suffix = text[last_idx]; // 접미사 추출

        uint8_t char_idx; // Occ 및 C 테이블에서 사용할 인덱스

        if (suffix == '$')
        {
            end_idx_onBWT = i; // 끝 표시 문자의 BWT에서의 인덱스 저장
            char_idx = 0; // $는 0번 인덱스
        }
        else
        {
            // A=1, C=2, G=3, T=4 (char_to_2bit 결과에 1을 더함)
            char_idx = char_to_2bit(suffix) + 1;
        }

        current_occ[char_idx]++;
        C_table[char_idx]++;

        // 현재 위치 i에 대한 Occ 값 기록
        for (uint8_t b = 0; b < 5; b++) Occ_table[b][i] = current_occ[b];

        // 4개씩 묶어서 1바이트로 압축 저장(2 * 4bit = 1byte)
        // 비트 연산으로 2비트씩 밀어넣음
        // (i % 4)가 0이면 첫 2비트, 1이면 다음 2비트...
        buffer |= (char_to_2bit(suffix) << (6 - (i % 4) * 2));

        if (i % 4 == 3 || i == length - 1)
        {
            // 완성된 1바이트를 벡터에 저장
            bwt_2bit.push_back(buffer); // ex. 10 11 11 01 -> 10111101 -> 189
            buffer = 0; // 버퍼 초기화
        }
    }


    // C-Table 완성 (개수 -> 시작 인덱스로 변환)
    size_t total = 0;
    for (uint8_t i = 0; i < 5; i++)
    {
        size_t count = C_table[i];
        C_table[i] = total;
        total += count;
    }

}
*/

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
        {
            positions.push_back(i);
        }
    }
    return positions;
}

/*
// LF Mapping 함수
// i번째 위치의 문자가 First Column에서 몇 번째인지 반환
size_t LF_mapping(size_t i, const vector<uint8_t>& bwt_2bit, size_t end_idx_onBWT, const vector<size_t>& C_table, const vector<vector<size_t>>& Occ_table)
{
    char c = get_bwt_char(bwt_2bit, i, end_idx_onBWT);

    uint8_t char_idx;
    if (c == '$') char_idx = 0;
    else char_idx = char_to_2bit(c) + 1;

    // LF(i) = C[c] + Occ[c][i]
    return C_table[char_idx] + Occ_table[char_idx][i];
}

// FM-Index 패턴 검색 함수
// 패턴이 reference에서 등장하는 모든 위치 반환
vector<size_t> FM_search(const string& pattern, const vector<uint8_t>& bwt_2bit, size_t end_idx_onBWT, const vector<size_t>& C_table, const vector<vector<size_t>>& Occ_table, const vector<size_t>& suffix_array)
{
    vector<size_t> positions; // 매핑된 위치들

    size_t length = suffix_array.size(); // BWT 길이 (원본 + $)

    // top, bot 포인터 초기화 (전체 범위)
    size_t top = 0;
    size_t bot = length;

    // 패턴을 역방향으로 검색
    for (int i = pattern.length() - 1; i >= 0; i--)
    {
        char c = pattern[i];

        uint8_t char_idx;
        if (c == '$') char_idx = 0;
        else char_idx = char_to_2bit(c) + 1;

        // top, bot 업데이트
        // Occ_table[char_idx][top-1] : top 이전까지의 누적 빈도
        size_t occ_top = (top > 0) ? Occ_table[char_idx][top - 1] : 0;
        size_t occ_bot = Occ_table[char_idx][bot - 1];

        top = C_table[char_idx] + occ_top;
        bot = C_table[char_idx] + occ_bot;

        // 매칭 없음
        if (top >= bot)
        {
            return positions; // 빈 벡터 반환
        }
    }

    // 매칭된 위치들 suffix array에서 찾기
    for (size_t i = top; i < bot; i++)
    {
        // suffix_array[i]는 $가 추가된 서열 기준이라 그대로 사용
        positions.push_back(suffix_array[i]);
    }

    return positions;
}
*/

int main() {
    // 파일에서 DNA 서열을 읽어옴
    string dna = load_fasta_of_len("chr22.fa", 10000); // chromosome 22
    //string dna = load_fasta("GCF_000001405.40_GRCh38.p14_genomic.fna"); // GRCh38 전체 유전체

    //cout << "추출된 DNA 길이: " << dna.length() << endl << endl;

    // 테스트용 짧은 DNA 서열
    //string dna = "AGCTACCAGGTG";
    string original_dna = dna; // 원본 서열 보존 (끝 표시 문자 추가 전)

    vector<uint64_t> bwt_2bit;

    /*
    *  끝 표시 문자의 인덱스 저장 변수
    *  2bit BWT에서는 끝 표시 문자를 00으로 처리하기 때문에, 별도의 인덱스 저장이 필요할 것 같아서 변수를 선언함
    *  불필요하다면 쓰지 않아도 됨
    */
    size_t end_idx_onBWT;

    vector<size_t> suffix_array;

    // $, A, C, G, T 순서
    // 각 문자의 BWT 내에서의 시작 인덱스 저장
    vector<size_t> C_table;
    // Occ_table[문자:0 ~ 4][인덱스:0 ~ length-1]에 해당하는 문자의 누적 빈도 저장
    // BWT 내에 각 인덱스에서의 A/C/G/T 누적 빈도 저장
    vector<vector<size_t>> Occ_table;

    // 2bit BWT 수행
    //BWT_2bit(dna, bwt_2bit, suffix_array, C_table, Occ_table, end_idx_onBWT);
    BWT_2bit(dna, bwt_2bit, suffix_array, C_table, Occ_table, end_idx_onBWT);

    /*
    // 2bit BWT에서 끝 표시 문자의 인덱스 출력
    cout << endl << endl;
    cout << "End marker index on BWT: " << end_idx_onBWT << endl;

    // SA 출력
    cout << "SA: ";
    for (size_t idx : suffix_array)
    {
        cout << idx << " ";
    }
    cout << endl << endl;

    // C-Table 출력
    for (size_t i = 0; i < C_table.size(); i++)
    {
        char c = (i == 0) ? '$' : (i == 1) ? 'A' : (i == 2) ? 'C' : (i == 3) ? 'G' : 'T';
        cout << "C[" << c << "] = " << C_table[i] << endl;
    }
    cout << endl;

    // Occ-Table 출력
    for (size_t i = 0; i < Occ_table.size(); i++)
    {
        char c = (i == 0) ? '$' : (i == 1) ? 'A' : (i == 2) ? 'C' : (i == 3) ? 'G' : 'T';
        cout << "Occ[" << c << "]: ";
        for (size_t j = 0; j < Occ_table[i].size(); j++)
        {
            cout << Occ_table[i][j] << " ";
        }
        cout << endl;
    }
    */

    /*
    // 패턴 검색 테스트
    string pattern = "AGG"; // 찾고 싶은 패턴
    //vector<size_t> result = FM_search(pattern, bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
    vector<size_t> result = FM_search_Optimized(pattern, bwt_packed, end_idx_onBWT, C_table, Occ_table, suffix_array);

    cout << "패턴 \"" << pattern << "\" 검색 결과: ";
    if (result.empty())
    {
        cout << "매칭 없음" << endl;
    }
    else
    {
        for (size_t pos : result)
        {
            cout << pos << " ";
        }
        cout << endl;
    }
    */

    // 1. reads 생성
    int M = 6400;
    int L = 16;
    vector<string> reads;
    srand(time(0));
    for (int i = 0; i < M; i++)
    {
        int start = rand() % (original_dna.length() - L);
        reads.push_back(original_dna.substr(start, L));
    }
    /*
    cout << "생성된 reads:" << endl;
    for (int i = 0; i < reads.size(); i++)
    {
        cout << "read[" << i << "]: " << reads[i] << endl;
    }
    cout << endl;
    */

    //// 2. reads 매핑
    //vector<vector<int>> mapping_table(dna.length(), vector<int>(4, 0));
    //cout << "reads 매핑 결과:" << endl;
    //for (int i = 0; i < reads.size(); i++)
    //{
    //    vector<size_t> positions = FM_search(reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
    //    if (positions.empty())
    //    {
    //        cout << "read[" << i << "] \"" << reads[i] << "\": 매칭 없음" << endl;
    //        continue;
    //    }
    //    for (size_t pos : positions)
    //    {
    //        cout << "read[" << i << "] \"" << reads[i] << "\": 위치 " << pos << endl;
    //        for (int j = 0; j < reads[i].length(); j++)
    //        {
    //            if (pos + j >= dna.length()) break;
    //            mapping_table[pos + j][char_to_2bit(reads[i][j])]++;
    //        }
    //    }
    //}
    //cout << endl;

    // 벤치마크 비교

    // FM-Index 시간 측정
    auto start_fm = chrono::high_resolution_clock::now();

    vector<vector<int>> mapping_table(original_dna.length(), vector<int>(4, 0));
    for (int i = 0; i < reads.size(); i++)
    {
        vector<size_t> positions = FM_search(reads[i], /*bwt_2bit*/bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
        for (size_t pos : positions)
        {
            for (int j = 0; j < reads[i].length(); j++)
            {
                if (pos + j >= original_dna.length()) break;
                mapping_table[pos + j][char_to_2bit(reads[i][j])]++;
            }
        }
    }

    auto end_fm = chrono::high_resolution_clock::now();
    auto duration_fm = chrono::duration_cast<chrono::microseconds>(end_fm - start_fm);
    cout << "FM-Index 실행시간: " << duration_fm.count() << " microseconds" << endl;


    // Trivial Sliding 시간 측정
    auto start_trivial = chrono::high_resolution_clock::now();

    vector<vector<int>> mapping_table_trivial(original_dna.length(), vector<int>(4, 0));
    for (int i = 0; i < reads.size(); i++)
    {
        vector<size_t> positions = Trivial_search(reads[i], original_dna, 0);
        for (size_t pos : positions)
        {
            for (int j = 0; j < reads[i].length(); j++)
            {
                if (pos + j >= original_dna.length()) break;
                mapping_table_trivial[pos + j][char_to_2bit(reads[i][j])]++;
            }
        }
    }

    auto end_trivial = chrono::high_resolution_clock::now();
    auto duration_trivial = chrono::duration_cast<chrono::microseconds>(end_trivial - start_trivial);
    cout << "Trivial 실행시간: " << duration_trivial.count() << " microseconds" << endl;


    // 3. Consensus 복원
    string recovered = "";
    for (int i = 0; i < original_dna.length(); i++)
    {
        int total = mapping_table[i][0] + mapping_table[i][1] + mapping_table[i][2] + mapping_table[i][3];
        if (total == 0)
        {
            recovered += original_dna[i]; // 리드가 매핑되지 않은 곳은 원본 서열을 그대로 복사
            continue;
        }
        int max_idx = 0;
        for (int j = 1; j < 4; j++)
        {
            if (mapping_table[i][j] > mapping_table[i][max_idx])
            {
                max_idx = j;
            }
        }
        recovered += bit_to_char(max_idx);
    }
    //cout << "원본 서열:  " << original_dna << endl;
    //cout << "복원 서열:  " << recovered << endl;

    // 4. 정확도 측정
    int match = 0;
    for (int i = 0; i < original_dna.length(); i++)
    {
        if (original_dna[i] == recovered[i]) match++;
        else cout << "불일치 위치: " << i << " (원본: " << original_dna[i] << ", 복원: " << recovered[i] << ")" << endl;
    }
    double accuracy = (double)match / original_dna.length() * 100;
    cout << "정확도: " << accuracy << "%" << endl;

    return 0;
}