# Progetto-MCS

## Relazione Part1

La cartella `Part1` contiene un progetto sviluppato per risolvere sistemi lineari di grandi dimensioni utilizzando diversi metodi iterativi. Il progetto sfrutta la libreria Eigen per la gestione efficiente di matrici dense e sparse.

### Obiettivo
L'obiettivo principale è confrontare le prestazioni e l'accuratezza di vari metodi iterativi per la risoluzione di sistemi lineari, sia su matrici dense che sparse. Vengono utilizzati dati di test forniti in formato Matrix Market (`.mtx`).

### Struttura del codice
- **main.cpp**: Punto di ingresso del programma. Carica le matrici di test, crea le versioni dense e sparse, e invoca la classe `linear_resolver` per eseguire i diversi metodi di risoluzione.
- **linear_resolver.hpp**: Definisce la classe template `linear_resolver`, che implementa i metodi iterativi Jacobi, Gauss-Seidel, Gradiente e Gradiente Coniugato. Gestisce la creazione dei vettori e il ciclo di risoluzione.
- **lr_utils.hpp**: Contiene funzioni di utilità per operazioni su matrici e vettori, come l'estrazione della diagonale e la sostituzione in avanti, sia per matrici dense che sparse.
- **dati/**: Cartella con i file di test in formato `.mtx`.
- **eigen-3.4.0/**: Versione locale della libreria Eigen necessaria per la compilazione.

### Funzionalità principali
- Supporto a matrici sia dense che sparse.
- Implementazione di metodi iterativi classici (Jacobi, Gauss-Seidel, Gradiente, Gradiente Coniugato).
- Possibilità di testare i metodi su matrici di esempio e su dati reali caricati da file.
- Utilizzo di Eigen per ottimizzare le operazioni su matrici e vettori.

### Compilazione ed esecuzione
La compilazione avviene tramite il `Makefile` presente nella cartella `Part1`. È sufficiente eseguire:

```
make
```

Per eseguire il programma:

```
make run
```

### Note aggiuntive
- Il codice è pensato per essere facilmente estendibile ad altri metodi iterativi.
- La modularità delle utility consente di adattare facilmente le funzioni a diversi tipi di matrice.