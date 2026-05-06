# Metody Obliczeniowe – Zadanie 11c

Numeryczne rozwiązanie równania dyfuzji metodami jawnymi i niejawnymi (Crank-Nicolson), z porównaniem błędów względem rozwiązania analitycznego.

---

## Opis problemu

Projekt rozwiązuje jednowymiarowe **równanie dyfuzji** (równanie ciepła):

$$\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2}$$

z warunkiem początkowym funkcji skokowej Heaviside'a:

$$u(x, 0) = \begin{cases} 1 & x < 0 \\ 0.5 & x = 0 \\ 0 & x > 0 \end{cases}$$

oraz warunkami brzegowymi Dirichleta na krańcach domeny `[-a, a]`.

Rozwiązanie analityczne wyraża się przez funkcję komplementarną błędu:

$$u(x, t) = \frac{1}{2} \, \mathrm{erfc}\!\left(\frac{x}{2\sqrt{Dt}}\right)$$

---

## Parametry symulacji

| Parametr | Wartość | Opis |
|---|---|---|
| `D` | 1.0 | Współczynnik dyfuzji |
| `t_max` | 2.0 | Czas symulacji [s] |
| `a` | `6·√(D·t_max)` | Połowa rozmiaru domeny |
| `λ` (jawna) | 0.4 | Liczba Couranta (< 0.5 → stabilna) |
| `λ` (niejawna) | 1.0 | Liczba Couranta dla metod niejawnych |

---

## Metody numeryczne

### 1. Metoda jawna (FTCS)
Prosty schemat różnic skończonych do przodu w czasie i centralnych w przestrzeni. Bezwarunkowo stabilna tylko dla `λ = D·dt/dx² ≤ 0.5`. Przy naruszeniu tego warunku program zgłasza wyjątek.

### 2. Metoda niejawna – algorytm Thomasa (Crank-Nicolson)
Schemat Crank-Nicolson prowadzi do układu równań z macierzą trójdiagonalną, rozwiązaną efektywnie **algorytmem Thomasa** (eliminacja Gaussa dla macierzy trójdiagonalnych). Bezwarunkowo stabilna.

### 3. Metoda niejawna – dekompozycja LU (Crank-Nicolson)
Ten sam schemat Crank-Nicolson, lecz układ rozwiązywany przez **rozkład LU** gęstej macierzy. Rozkład wykonywany jednorazowo przed pętlą czasową – każdy krok czasowy wymaga tylko podstawiania w przód i wstecz. Bezwarunkowo stabilna.

---

## Struktura plików

```
.
├── main.cpp          # Punkt wejścia; analiza błędów w czasie
├── params.h          # Globalne stałe fizyczne i numeryczne
├── analytical.cpp/h  # Rozwiązanie analityczne (erfc)
├── CALERF.cpp/h      # Implementacja funkcji erfc wysokiej precyzji (long double)
├── explicit.cpp/h    # Metoda jawna (FTCS)
├── implicit.cpp/h    # Metody niejawne: Thomas i LU
├── linalg.cpp/h      # Algorytm Thomasa + rozkład LU + podstawianie
├── output.cpp/h      # Zapis wyników do plików .dat
└── CMakeLists.txt    # Konfiguracja budowania (CMake, C++20)
```

---

## Budowanie i uruchamianie

### Wymagania
- Kompilator C++20 (GCC ≥ 10, Clang ≥ 12, MSVC ≥ 2019)
- CMake ≥ 3.31

### Kompilacja

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

### Uruchomienie

```bash
./untitled2
```

Program przeprowadza analizę błędów dla trzech metod i zapisuje wyniki do plików:

```
time_error_explicit.dat
time_error_thomas.dat
time_error_lu.dat
```

Każdy plik zawiera kolumny `t` i `max_error` – maksymalny błąd bezwzględny względem rozwiązania analitycznego w danej chwili czasu.

---

## Wyjście programu

```
[e] explicit...
t = 0.1, max error = ...
...

[t] thomas...
t = 0.1, max error = ...
...

[l] lu...
t = 0.1, max error = ...
...

Analiza zakończona. Dane zapisane do plików:
 - time_error_explicit.dat
 - time_error_thomas.dat
 - time_error_lu.dat

Total execution time: N seconds
```

---

## Uwagi implementacyjne

- Funkcja `erfc` obliczana jest z precyzją `long double` przy użyciu biblioteki **CALERF** (Cody, 1969), co zapewnia dokładność rozwiązania analitycznego na poziomie maszynowym.
- Rozkład LU macierzy trójdiagonalnej jest wykonywany **poza pętlą czasową**, co eliminuje wielokrotne faktoryzacje i przyspiesza obliczenia kosztem większego zużycia pamięci względem algorytmu Thomasa.
- Warunek stabilności metody jawnej (`λ ≤ 0.5`) jest sprawdzany w czasie wykonania – przekroczenie rzuca `std::invalid_argument`.
