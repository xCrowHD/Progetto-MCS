# Progetto 1: Mini libreria per sistemi lineari (Parte I)

## 1. Introduzione

In questo documento presentiamo la **Parte I** del progetto di realizzazione di una mini libreria per la risoluzione di sistemi lineari mediante metodi iterativi. L’obiettivo è implementare e validare i **quattro** solutori iterativi per matrici simmetriche e definite positive:

1. Metodo di Jacobi
2. Metodo di Gauss–Seidel
3. Metodo del Gradiente
4. Metodo del Gradiente Coniugato

---

## 2. Nozioni teoriche

### 2.1 Matrici simmetriche e definite positive

Le metodologie iterative qui descritte richiedono matrici $A \in \mathbb{R}^{n\times n}$ tali che:

* **Simmetria**: $A = A^T$
* **Definità positiva**: $x^T A x > 0$ per ogni $x \neq 0$

### 2.2 Metodo di Jacobi

**Descrizione**: Il metodo di Jacobi separa la matrice $A$ in una parte diagonale $D$ e nel resto $R = A - D$. A partire da un’approssimazione iniziale $x^{(0)}$, la regola iterativa è:

$$
x^{(k+1)} = D^{-1} \bigl(b - R\,x^{(k)}\bigr)
$$

La convergenza è garantita se $A$ è simmetrica e definita positiva.

### 2.3 Metodo di Gauss–Seidel

**Descrizione**: Il metodo di Gauss–Seidel sfrutta i valori aggiornati di $x^{(k+1)}$ durante la stessa iterazione. Si esprime $A = D + L + U$ (triangolare inferiore $L$, superiore $U$) e si applica:

$$
x^{(k+1)} = (D + L)^{-1} \bigl(b - U\,x^{(k)}\bigr)
$$

Rispetto a Jacobi, converge generalmente più velocemente per matrici SPD.

### 2.4 Metodo del Gradiente

**Descrizione**: Basato sulla minimizzazione della funzione quadratica $f(x) = \frac{1}{2} x^T A x - b^T x$. Definito il residuo $r^{(k)} = b - A x^{(k)}$, l’iterazione segue:

$$
x^{(k+1)} = x^{(k)} + \alpha^{(k)} r^{(k)},
$$

dove

$$
\alpha^{(k)} = \frac{r^{(k)T} r^{(k)}}{r^{(k)T} A r^{(k)}}.
$$

### 2.5 Metodo del Gradiente Coniugato

**Descrizione**: Costruisce direzioni $p^{(k)}$ coniugate rispetto a $A$. Algoritmo:

1. $x^{(0)}$ iniziale, $r^{(0)} = b - A x^{(0)}$, $p^{(0)} = r^{(0)}$.
2. Per $k=0,1,\dots$ finché non converge:

   * $\alpha^{(k)} = \frac{r^{(k)T} r^{(k)}}{p^{(k)T} A p^{(k)}}$
   * $x^{(k+1)} = x^{(k)} + \alpha^{(k)} p^{(k)}$
   * $r^{(k+1)} = r^{(k)} - \alpha^{(k)} A p^{(k)}$
   * $\beta^{(k)} = \frac{r^{(k+1)T} r^{(k+1)}}{r^{(k)T} r^{(k)}}$
   * $p^{(k+1)} = r^{(k+1)} + \beta^{(k)} p^{(k)}$

Questo metodo converge in al massimo $n$ iterazioni per una matrice $n\times n$ SPD.

---

## 3. Struttura del codice

La cartella principale di questa parte del progetto è `Part1`, organizzata come segue:

```
Part1/
├── dati/                   # Contiene matrici in formato .mtx e .csv
├── linear_resolver.hpp     # Implementazione dei metodi iterativi (Jacobi, Gauss–Seidel, Gradiente, Gradiente Coniugato)
├── lr_utils.hpp            # Funzioni di supporto per LinearResolver
├── lr_test.hpp             # Generazione di tabelle a terminale con le seguenti colonne:
│   │ Nome Metodo │ Numero iterazioni │ Residuo relativo │ Errore relativo │ Errore assoluto │ Tempo (s)
├── main.cpp                # Interfaccia a terminale per:
│   │ • Eseguire tutti i test su file .mtx con varie tolleranze
│   │ • Selezionare una singola matrice (.mtx o .csv) e una tolleranza per visualizzare una tabella
├── Makefile                # Target:
│   │ • `make build` : compila il progetto
│   │ • `make run`   : esegue l’eseguibile con `main.cpp`
├── tabulate.hpp            # Libreria header-only per la formattazione delle tabelle

# Directory esterna
├── ../Eigen/               # Libreria Eigen per strutture dati dense e sparse (inclusa come dipendenza locale)
```

## 4. Classi e componenti principali

### 4.1 Classe `linear_resolver<MatrixType>`

Questa è la classe principale che implementa i **metodi risolutivi iterativi** per sistemi lineari della forma `Ax = b`. È una classe template che può lavorare sia con matrici dense (`Eigen::MatrixXd`) che sparse (`Eigen::SparseMatrix<double>`), grazie alla flessibilità della libreria Eigen.

#### **Costruttore**

```cpp
linear_resolver(const MatrixType& A, const double& tol)
```

* Inizializza il sistema con una matrice `A` e una tolleranza `tol`.
* Viene generato un vettore soluzione noto `x` come vettore di `1`, e il vettore dei termini noti `b` viene calcolato come `A * x`.

#### **Metodi principali**

La classe implementa quattro algoritmi di risoluzione iterativa, ognuno dei quali restituisce una tupla con:

```cpp
(x, numero_iterazioni, residuo_relativo, tempo_sec, errore_assoluto, errore_relativo)
```

* `jacobi_resolver(...)`
  Implementa il metodo di Jacobi classico.
  Utilizza l'inverso della diagonale `D⁻¹` della matrice per aggiornare la soluzione iterativamente.

* `gauss_resolver(...)`
  Implementa il metodo di **Gauss-Seidel** tramite una soluzione approssimata con **forward substitution** sulla matrice triangolare inferiore.

* `gradient_resolver(...)`
  Risolve il sistema tramite **metodo del gradiente** (discesa del gradiente), calcolando la direzione del massimo calo.

* `conjugate_gradient_resolver(...)`
  Implementazione completa del **metodo del gradiente coniugato**, particolarmente efficiente su sistemi simmetrici definiti positivi.

#### **Metodi di accesso**

* `getA() const` — Ritorna la matrice `A`.
* `getb() const` — Ritorna il vettore `b`.

#### **Dettagli aggiuntivi**

* I metodi misurano automaticamente il tempo di esecuzione usando `std::chrono`.
* Il numero massimo di iterazioni è impostato a 20000.
* La classe è completamente generica e modulare, pensata per essere estendibile e riutilizzabile.

---

### 4.2 Classe `lr_utils`

Questo namespace raccoglie una serie di funzioni di utilità per supportare l’implementazione dei metodi iterativi nella classe `linear_resolver`. Le principali funzionalità includono:

* **Estrazione della diagonale**:
  Funzioni `getD` sovraccaricate per estrarre il vettore diagonale principale da matrici dense (`Eigen::MatrixXd`) o sparse (`Eigen::SparseMatrix<double>`).

* **Sostituzione in avanti (Forward substitution)**:
  Implementazioni personalizzate di forward substitution per risolvere sistemi triangolari inferiori $P y = r$, adattate sia per matrici dense sia sparse.
  Questo perché le funzioni predefinite di Eigen di forward substituition non sono ottimizzate per il nostro use case.

* **Estrazione della matrice triangolare inferiore**:
  Funzioni `getLowerTriangular` per ottenere la matrice triangolare inferiore (inclusa la diagonale) da una matrice densa o sparsa.
  Questa matrice $P = D + L$ viene usata nel metodo di Gauss–Seidel.

**Dettagli implementativi importanti**:

* La forward substitution per matrici sparse viene eseguita convertendo temporaneamente la matrice in formato `RowMajor` per iterare efficientemente riga per riga.
* Per la matrice triangolare inferiore sparsa, viene creata tramite un vettore di triplette (coordinate non nulle) per mantenere la struttura sparsa e ottimizzare le operazioni successive.

---

### 4.3 Classe `lr_test`

Classe template `lr_test<MatrixType>` che si occupa di testare i metodi iterativi implementati in `linear_resolver` su matrici caricate da file. Supporta sia matrici dense (`Eigen::MatrixXd`) che sparse (`Eigen::SparseMatrix<double>`).

* **Caricamento matrici automatico da cartella**:
  La classe esegue una scansione automatica della cartella `./dati` per individuare tutti i file con estensione `.mtx`.
  Questo significa che non ci sono matrici "hardcoded" nel codice:

  * Per aggiungere nuovi test è sufficiente inserire nuovi file `.mtx` nella cartella `dati`
  * Per rimuovere test basta cancellare o spostare i file `.mtx`
    Non è necessario modificare il codice per cambiare i dataset testati, garantendo flessibilità e facilità di estensione.


* **Test multipli su tolleranze**:
  Per ogni matrice e per un insieme di tolleranze predefinite (`1e-4`, `1e-6`, `1e-8`, `1e-10`), applica tutti e quattro i metodi iterativi.

* **Creazione tabelle di risultati**:
  Per ogni esecuzione genera tabelle formattate con la libreria `tabulate` contenenti:

  * Nome metodo
  * Numero di iterazioni
  * Residuo relativo $\|b - Ax\| / \|b\|$
  * Errore relativo $\|x - x_{\text{exact}}\| / \|x_{\text{exact}}\|$
  * Errore assoluto $\|x - x_{\text{exact}}\|$
  * Tempo di esecuzione (secondi)

* **Colorazione automatica**:
  Evidenzia in verde i valori migliori (minimi) in ciascuna colonna (es. minimo numero di iterazioni, minimo residuo, ecc).

* **Statistiche aggregate**:
  Dopo aver eseguito tutti i test, crea una tabella di statistiche che conta quante volte ciascun metodo ha ottenuto il valore minimo per ogni metrica.

**Dettagli implementativi importanti**:

* Usa C++17, `std::filesystem` per navigare i file.
* Usa la libreria `Eigen` per manipolare matrici sparse e dense.
* Implementa metodi per estrarre i dati dai risultati e formattarli per la visualizzazione.
* La gestione del colore e dello stile è gestita tramite `tabulate`.

---

### 4.4 Classe `main.cpp`

la classe `main.cpp` implementa un'interfaccia a terminale interattiva per eseguire i test della mini libreria.

**Funzionalità principali:**

* **Menu semplice** con tre opzioni:

  1. Eseguire automaticamente i test su tutte le matrici sparse `.mtx` presenti nella cartella `./dati`.
  2. Eseguire il test su una singola matrice scelta dall’utente, fornendo il percorso di un file `.mtx` (sparsa) o `.csv` (densa).
  3. Uscire dal programma.

* **Gestione dinamica dei file**:

  * Nella prima opzione, il programma scansiona automaticamente la directory `./dati` per trovare tutti i file `.mtx` senza hardcoding, permettendo di aggiungere o rimuovere matrici senza modifiche al codice.
  * Nella seconda opzione, è possibile specificare manualmente un file di matrice da testare, sia in formato `.mtx` che `.csv`.

* **Caricamento matrici**:

  * Le matrici sparse `.mtx` sono caricate con `Eigen::loadMarket`.
  * Le matrici dense `.csv` sono lette tramite parsing personalizzato (supportano separatore specificato).

* **Input utente**:

  * L’utente può inserire la tolleranza desiderata per la risoluzione iterativa, con validazione del valore.

* **Output**:

  * I risultati sono mostrati in tabelle formattate con la libreria `tabulate.hpp`, evidenziando statistiche come numero di iterazioni, residuo relativo, errori e tempi di esecuzione per ciascun metodo.

---
