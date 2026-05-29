#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// 모든 ASCII 문자를 커버하기 위한 알파벳 크기 고정 (256개)
const int ALPHABET_SIZE = 256;

// 대용량 데이터 처리 시 메모리 폭발을 방지하기 위한 샘플링 간격 설정
// 모든 인덱스 대신 64칸마다 한 번씩만 Occ 값을 기록하여 메모리를 1/64로 절감
const int CHECKPOINT = 64;

// =========================================================================
// [SA-IS 하위 시스템] 접미사 유도 정렬 (Induced Sorting)
// =========================================================================
void induce(const vector<int>& s, vector<int>& sa, const vector<bool>& ls,
            const vector<int>& bucket_start, const vector<int>& bucket_end) {
    int n = s.size();
    
    // 1. L-type 접미사 유도 (정방향 스캔: 앞에서부터 뒤로)
    vector<int> buf = bucket_start;
    for (int i = 0; i < n; i++) {
        int v = sa[i] - 1;
        // 내 앞의 접미사가 존재하고, 그것이 L-type인 경우 버킷의 앞에서부터 채움
        if (v >= 0 && !ls[v]) {
            sa[buf[s[v]]++] = v;
        }
    }
    
    // 2. S-type 접미사 유도 (역방향 스캔: 뒤에서부터 앞으로)
    buf = bucket_end;
    for (int i = n - 1; i >= 0; i--) {
        int v = sa[i] - 1;
        // 내 앞의 접미사가 존재하고, 그것이 S-type인 경우 버킷의 뒤에서부터 채움
        if (v >= 0 && ls[v]) {
            sa[--buf[s[v]]] = v;
        }
    }
}

// =========================================================================
// [SA-IS 핵심 모듈] O(N) 선형 시간 접미사 배열 생성
// =========================================================================
vector<int> sa_is(const vector<int>& s, int upper) {
    int n = s.size();
    if (n == 0) return {};
    if (n == 1) return {0};
    
    vector<int> sa(n, -1);
    vector<bool> ls(n);
    
    // 문자열 끝에서부터 거꾸로 읽으며 L-type(대) / S-type(소) 분류
    ls[n - 1] = true; // Sentinel 문자(종단 기호)는 사전순 최솟값이므로 항상 S-type
    for (int i = n - 2; i >= 0; i--) {
        ls[i] = (s[i] < s[i + 1] || (s[i] == s[i + 1] && ls[i + 1]));
    }

    // 각 문자별 빈도수 측정 및 버킷 범위(시작점/끝점) 계산
    vector<int> bucket_cnt(upper + 1, 0);
    for (int x : s) bucket_cnt[x]++;
    vector<int> bucket_start(upper + 1, 0), bucket_end(upper + 1, 0);
    int sum = 0;
    for (int i = 0; i <= upper; i++) {
        bucket_start[i] = sum;
        sum += bucket_cnt[i];
        bucket_end[i] = sum;
    }

    // LMS(Leftmost S-type) 위치 탐색 (L-type 다음에 등장하는 최초의 S-type 마킹)
    vector<int> lms_pos;
    vector<int> lms_map(n, -1);
    int m = 0;
    for (int i = 1; i < n; i++) {
        if (!ls[i - 1] && ls[i]) {
            lms_map[i] = m++;
            lms_pos.push_back(i);
        }
    }

    // [1차 추측 정렬] LMS를 대략적인 버킷 끝자리에 배치 후 유도 정렬 수행
    vector<int> buf = bucket_end;
    for (int p : lms_pos) {
        sa[--buf[s[p]]] = p;
    }
    induce(s, sa, ls, bucket_start, bucket_end);

    // LMS 블록들의 실제 정렬 순서를 확인하여 압축 문자열 S1 생성 절차 시작
    if (m > 0) {
        vector<int> sorted_lms;
        for (int x : sa) {
            if (lms_map[x] != -1) sorted_lms.push_back(x);
        }

        int num_names = 0;
        vector<int> lms_names(m, -1);
        lms_names[lms_map[sorted_lms[0]]] = 0;

        // 인접한 LMS 블록끼리 세부 비교하여 고유 번호(Name) 부여
        for (int i = 1; i < m; i++) {
            int u = sorted_lms[i - 1], v = sorted_lms[i];
            bool same = true;
            for (int d = 0; ; d++) {
                if (s[u + d] != s[v + d] || ls[u + d] != ls[v + d]) {
                    same = false;
                    break;
                }
                if (d > 0 && ((!ls[u + d - 1] && ls[u + d]) || (!ls[v + d - 1] && ls[v + d]))) {
                    break; // 다음 LMS 캐릭터 경계면에 닿으면 비교 종료
                }
            }
            if (!same) num_names++;
            lms_names[lms_map[v]] = num_names;
        }

        // 압축된 새로운 정수형 문자열 집합 s1 구축
        vector<int> s1(m);
        for (int i = 0; i < m; i++) {
            s1[i] = lms_names[lms_map[lms_pos[i]]];
        }

        // 재귀 호출 분기: 고유 번호가 부족하면 완전히 정렬될 때까지 재귀 호출
        vector<int> sa1;
        if (num_names + 1 == m) {
            sa1.resize(m);
            for (int i = 0; i < m; i++) sa1[s1[i]] = i;
        } else {
            sa1 = sa_is(s1, num_names);
        }

        // [2차 최종 정렬] 재귀 결과로 얻은 정확한 LMS 순서를 적용하여 완벽한 유도 정렬 수행
        fill(sa.begin(), sa.end(), -1);
        buf = bucket_end;
        for (int i = m - 1; i >= 0; i--) {
            int p = lms_pos[sa1[i]];
            sa[--buf[s[p]]] = p;
        }
        induce(s, sa, ls, bucket_start, bucket_end);
    }
    return sa;
}

// SA-IS 알고리즘 호출을 편하게 도와주는 래퍼 인터페이스
vector<int> build_suffix_array(const string& text) {
    vector<int> s;
    int upper = 0;
    // 0번 값을 Sentinel로 쓰기 위해 일반 문자들의 값을 1씩 증가시켜 인코딩
    for (char c : text) {
        s.push_back((unsigned char)c + 1);
        upper = max(upper, s.back());
    }
    s.push_back(0); // 문자열 맨 끝에 최솟값 0(종단 문자 '$') 추가
    return sa_is(s, upper);
}

// =========================================================================
// [핵심 최적화] 단일 패스(Single Pass) BWT 생성 및 FM-Index 테이블 일괄 구축
// =========================================================================
void build_fm_index_single_pass(const string& text, const vector<int>& sa, 
                                string& bwt, vector<int>& C, vector<vector<int>>& Occ_sampled) {
    // 원본 문자열에 종단 기호를 붙여 실제 연산 데이터 크기 획정
    string text_with_sentinel = text + "$";
    int n = text_with_sentinel.size();

    // ---------------------------------------------------------------------
    // STEP 1: C 테이블(Count Table) 빌드
    // ---------------------------------------------------------------------
    // 전체 문자열을 딱 한 번만 순회하여 전체 빈도수를 확보 ($0,000,000 기준 약 0.02초)
    vector<int> counts(ALPHABET_SIZE, 0);
    for (char c : text_with_sentinel) {
        counts[(unsigned char)c]++;
    }
    
    // 사전순으로 나보다 작은 문자들의 누적 빈도수를 계산하여 C 배열 완성
    int sum = 0;
    C.assign(ALPHABET_SIZE, 0);
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        C[i] = sum;
        sum += counts[i];
    }

    // ---------------------------------------------------------------------
    // STEP 2: 단일 패스 BWT 추출 및 Occ 테이블 샘플링 기록
    // ---------------------------------------------------------------------
    bwt.resize(n);
    int sampled_rows = (n / CHECKPOINT) + 1; // 64칸 간격으로 저장할 행 크기 계산
    Occ_sampled.assign(ALPHABET_SIZE, vector<int>(sampled_rows, 0));

    // 루프를 돌며 실시간으로 각 문자의 등장 횟수를 누적해 나갈 변수
    vector<int> current_counts(ALPHABET_SIZE, 0);

    // 접미사 배열(SA)을 순회하면서 BWT를 뽑는 동시에 카운팅 수행
    for (int i = 0; i < n; i++) {
        // 정렬된 접미사 시작 인덱스의 '바로 직전 문자'가 BWT 문자가 됨
        int idx = (sa[i] - 1 + n) % n;
        char c = text_with_sentinel[idx];
        bwt[i] = c;

        // BWT 문자열이 쌓이는 순서대로 실시간 카운트 가산
        current_counts[(unsigned char)c]++;

        // 대용량 처리를 위해 매 인덱스마다 기록하지 않고, i+1이 64의 배수가 될 때만 스냅샷 저장
        // 이로 인해 내부에 매번 256회씩 돌던 낭비 루프가 제거됨
        if ((i + 1) % CHECKPOINT == 0) {
            int chk_idx = (i + 1) / CHECKPOINT;
            for (int alpha = 0; alpha < ALPHABET_SIZE; alpha++) {
                Occ_sampled[alpha][chk_idx] = current_counts[alpha];
            }
        }
    }
}

// =========================================================================
// [조회 모듈] 띄엄띄엄 저장된 Occ 테이블에서 원하는 범위의 카운트를 계산하는 함수
// =========================================================================
int get_occ(const string& bwt, const vector<vector<int>>& Occ_sampled, char c, int i) {
    if (i <= 0) return 0;

    // 1. 요청한 인덱스 i와 가장 가까운 전방의 체크포인트 블록 위치 계산
    int chk_idx = i / CHECKPOINT;
    int start_pos = chk_idx * CHECKPOINT;
    
    // 2. 이미 기록해 둔 체크포인트 시점의 카운트 값을 베이스로 가져옴 (O(1))
    int count = Occ_sampled[(unsigned char)c][chk_idx];

    // 3. 체크포인트 지점부터 실제 타겟 인덱스 i 직전까지 남은 자투리(최대 63칸)만 부분 스캔
    // 대용량 데이터 환경에서도 최대 연산 횟수가 63회로 제한되므로 상수의 속도 보장
    for (int k = start_pos; k < i; k++) {
        if (bwt[k] == c) {
            count++;
        }
    }
    return count;
}

// =========================================================================
// 메인 실행 흐름
// =========================================================================
int main() {
    string text = "banana";
    
    // 1단계: 선형 시간 내에 Suffix Array 확보
    vector<int> sa = build_suffix_array(text);
    
    // 2단계: 추가적인 다중 순회 비용 없이, 한 바퀴만 돌아서 BWT와 FM-Index용 테이블을 일괄 계산
    string bwt;
    vector<int> C;
    vector<vector<int>> Occ_sampled;
    build_fm_index_single_pass(text, sa, bwt, C, Occ_sampled);

    // 인덱스 구축 결과 및 샘플링 쿼리 정상 작동 검증
    cout << "BWT Result : " << bwt << "\n";
    cout << "C['a']     : " << C['a'] << "\n";
    
    // FM-Index 기반 문자열 검색 시 역방향 탐색 공식에 들어가는 핵심 쿼리 기능 수행 테스트
    // 0번째부터 4번째(즉, 5개 문자 범위)까지 문자 'a'가 몇 번 나왔는지 정확히 판별
    cout << "Occ('a', 5): " << get_occ(bwt, Occ_sampled, 'a', 5) << "\n";

    return 0;
}