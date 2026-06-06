#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <bitset>
#include <cstdlib> 
#include <ctime>  
#include <chrono>


using namespace std;

string load_fasta(const string& filename)
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

	while (getline(file, line)) // 파일에서 한 줄씩 읽기
    {
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
            }
            // 'N'은 무시하고 저장하지 않음
        }
    }

    file.close();

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
        case '$': return 0; // 끝 표시 문자는 00으로 처리 (특수 인덱스로 구분)
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


// 2bit BWT에서 i번째 문자를 꺼내는 함수
char get_bwt_char(const vector<uint8_t>& bwt_2bit, size_t i, size_t end_idx)
{
    if (i == end_idx) return '$';

    size_t byte_idx = i / 4;
    size_t bit_offset = i % 4;

    uint8_t bits = (bwt_2bit[byte_idx] >> (6 - bit_offset * 2)) & 0x03;

    return bit_to_char(bits);
}

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

/*
// 2비트 BWT를 파일로 저장하는 함수 (메타데이터로 원본 서열의 길이도 함께 저장)
void save_bwt_2bit(const string& bwt, const string& filename)
{
    std::ofstream fout(filename, std::ios::binary);

    // 1. 원본 서열의 길이를 먼저 저장 (복구할 때 필요)
    uint64_t length = bwt.length();
    fout.write(reinterpret_cast<const char*>(&length), sizeof(length)); // 메타데이터로 길이 저장

	// 2. 4개씩 묶어서 1바이트로 압축 저장 (2 * 4bit =  1byte)
    unsigned char buffer = 0;
    for (size_t i = 0; i < length; ++i)
    {
        // 비트 연산으로 2비트씩 밀어넣음
        // (i % 4)가 0이면 첫 2비트, 1이면 다음 2비트...
        buffer |= (char_to_2bit(bwt[i]) << (6 - (i % 4) * 2));

        if (i % 4 == 3 || i == length - 1)
        {
            fout.put(buffer);
            buffer = 0; // 버퍼 초기화
        }
    }
    fout.close();
}
*/


// Trivial Sliding 함수 (벤치마크용)
vector<size_t> Trivial_search(const string& read, const string& reference, int max_mismatch)
{
    vector<size_t> positions;
    int L = read.length();
    int N = reference.length();

    for (int i = 0; i <= N - L; i++)
    {
        int mismatch = 0;
        for (int j = 0; j < L; j++)
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




int main() {
	// 파일에서 DNA 서열을 읽어옴
    //string dna = load_fasta("chr22.fa"); // chromosome 22
	//string dna = load_fasta("GCF_000001405.40_GRCh38.p14_genomic.fna"); // GRCh38 전체 유전체

    //cout << "추출된 DNA 길이: " << dna.length() << endl << endl;

    /*
	int print_length = 10000; // 최대 10000자까지만 출력
    for(char c : dna)
    {
        cout << c;
        if (--print_length == 0)
        {
            cout << "\n... (생략) ...\n";
            break;
		}
	}
    */
    
	// 테스트용 짧은 DNA 서열
	string dna = "AGCTACCAGGTG";
	string original_dna = dna; // 원본 서열 보존 (끝 표시 문자 추가 전)

	string bwt;
	vector<uint8_t> bwt_2bit;

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

    /*
    // BWT 수행
	BWT(dna, bwt, suffix_array, C_table, Occ_table);
	// BWT 결과 출력
    cout << "BWT: " << bwt << endl << endl;
    */
    
	// 2bit BWT 수행
	BWT_2bit(dna, bwt_2bit, suffix_array, C_table, Occ_table, end_idx_onBWT);

	// 2bit BWT 결과 출력
	cout << "BWT 2bit: ";
	for (uint8_t byte : bwt_2bit) cout << bitset<8>(byte) << " "; // 2bit로 압축된 BWT를 8비트 이진수로 출력 (각 바이트는 4개의 염기를 나타냄)

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


    // 패턴 검색 테스트
    string pattern = "AGG"; // 찾고 싶은 패턴
    vector<size_t> result = FM_search(pattern, bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);

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

    

    // 1. reads 생성
    int M = 10;
    int L = 4;
    vector<string> reads;
    srand(time(0));
    for (int i = 0; i < M; i++)
    {
        int start = rand() % (original_dna.length() - L);
        reads.push_back(original_dna.substr(start, L));
    }
    cout << "생성된 reads:" << endl;
    for (int i = 0; i < reads.size(); i++)
    {
        cout << "read[" << i << "]: " << reads[i] << endl;
    }
    cout << endl;

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
        vector<size_t> positions = FM_search(reads[i], bwt_2bit, end_idx_onBWT, C_table, Occ_table, suffix_array);
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
            recovered += original_dna[i];
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
    cout << "원본 서열:  " << original_dna << endl;
    cout << "복원 서열:  " << recovered << endl;

    // 4. 정확도 측정
    int match = 0;
    for (int i = 0; i < original_dna.length(); i++)
    {
        if (original_dna[i] == recovered[i]) match++;
    }
    double accuracy = (double)match / original_dna.length() * 100;
    cout << "정확도: " << accuracy << "%" << endl;

    return 0;  
}