# GitHub 업로드 가이드 (처음 사용자용)

## 1단계: GitHub 계정 만들기 (5분)

1. 브라우저에서 **https://github.com** 접속
2. 우측 상단 **Sign up** 클릭
3. 이메일, 비밀번호, 사용자명 입력
   - 사용자명 추천: `sungsoo-jung` 또는 `sjung-rheum`
4. 이메일 인증 완료

---

## 2단계: 새 Repository 만들기 (2분)

1. 로그인 후 우측 상단 **+** 버튼 → **New repository**
2. 설정:
   - **Repository name**: `3AIF-ODE-model`
   - **Description**: `ODE simulation code for the 3-Axis Integrative Framework (3-AIF)`
   - **Public** 선택 (논문 공개용)
   - ✅ Add a README file — **체크하지 마세요** (이미 만들어 놓음)
   - .gitignore, License — **None으로 두세요** (이미 만들어 놓음)
3. **Create repository** 클릭

---

## 3단계: 파일 업로드 (3분)

빈 repository 화면이 나타나면:

1. **uploading an existing file** 링크 클릭 (또는 **Add file** → **Upload files**)
2. 다운로드 받은 파일들을 **드래그 앤 드롭**:
   - `3aif_simulation.py`
   - `README.md`
   - `LICENSE`
   - `requirements.txt`
   - `.gitignore`
3. 하단 **Commit changes** 영역:
   - 메시지: `Initial commit: 3-AIF ODE simulation code`
4. **Commit changes** 버튼 클릭

---

## 4단계: 원고 Data Availability 수정

현재:
> ODE model code and parameters are available from the corresponding author on request.

변경:
> ODE model code, default parameter values, and simulation scripts are publicly available at https://github.com/[사용자명]/3AIF-ODE-model

---

## 5단계 (선택): Zenodo DOI 발급

영구 DOI가 필요하면:

1. **https://zenodo.org** 접속 → **Log in with GitHub**
2. 우측 상단 드롭다운 → **GitHub**
3. `3AIF-ODE-model` repository 옆 토글 **ON**
4. GitHub로 돌아가서 repository → **Releases** → **Create a new release**
   - Tag: `v1.0`
   - Title: `3-AIF ODE Model v1.0 (bioRxiv submission)`
5. **Publish release** 클릭
6. Zenodo가 자동으로 DOI 생성 (몇 분 소요)
7. 생성된 DOI를 원고에 추가

---

## 문제 해결

| 문제 | 해결 |
|------|------|
| `.gitignore`가 안 보임 | 파일 이름이 점(.)으로 시작하면 Mac/Windows에서 숨겨짐. 터미널에서 확인하거나 GitHub 웹에서 직접 생성 |
| 파일이 너무 큼 | GitHub은 파일당 100MB 제한. 시뮬레이션 코드는 문제없음 |
| 수정이 필요함 | GitHub 웹에서 파일 클릭 → 연필(✏️) 아이콘 → 수정 → Commit |

---

## Repository 구조 요약

```
3AIF-ODE-model/
├── README.md              ← 프로젝트 설명 (영어)
├── LICENSE                ← MIT 라이선스
├── requirements.txt       ← Python 패키지 목록
├── .gitignore             ← Git 제외 파일 목록
└── 3aif_simulation.py     ← 핵심 시뮬레이션 코드
```

실행 방법:
```bash
pip install -r requirements.txt
python 3aif_simulation.py
```

→ `figures/` 폴더에 Figure 2, 3B, 3C, 4A 이미지 생성
