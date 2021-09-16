# AntiCov
AxCalc.py jest skryptem pythona pozwalającym na wygenerowanie całego zbioru inputów do Virtual AFM. Docelowo będzie uruchamiany automatycznie przez przeglądarkę, ale na razie trzeba go dokładnie przetestować.

Przed uruchomieniem
--------------------------
Skrypt napisany jest w standardzie Python3 (ogólnie obowiązujący). Wymaga jednak instalacji kilku dodatkowych pakietów pythona, nie zawartych w standardowej instalacji. Są to NumPy, SciPy, oraz MDAnalysis.
Instalujemy je używając komendy pip:

 _pip install numpy scipy MDAnalysis_

lub za pomocą conda:

_conda install numpy scipy MDAnalysis_

Uruchamianie
------------------------
program uruchamiamy w terminalu (linux) wpisując w wierszu poleceń (będąc w katalogu w którym ściągnęliśmy skrypt):

_python AxCalc.py_  (lub python3 AxCalc.py  jeżeli python3 nie jest domyślną instancją pythona.) i w tej samej lini podajemy parametry wejściowe.

Program wywołujemy podając mu następujące parametry:

_python AxCalc.py file.pdb file.psf file.vel file.coor file.xsc toppar.zip template.inp template.run ‘selection constraints’ ‘selection pull’_

file.pdb        - struktura pdb naszego układu

file.psf        - struktura psf układu

file.vel        - informacja o prędkościach (zakładamy, że układ jest już po ekwilibracji)  

file.coor    - informacja o aktualnych współrzędnych wszystkich atomów (j.w)

file.xsc        - contains the periodic cell parameters and extended system variables

toppar.zip    - jeżeli mamy tylko jeden plik z parametrami Charmm to wpisujemy go tutaj
  (np. param1.inp). Jeżeli jednak tych plików z parametrami jest więcej, trzeba 
  włożyć je do folderu toppar i “zzipować” (zip -r toppar.zip toppar)

template.inp    - tutaj mamy plik wejściowy do namd, w którym ustawiamy sobie wszystkie 
  potrzebne nam parametry symulacji. Na podstawie tego pliku będą generowane poszczególne inputy do SMD więc ważne jest, by część opisująca SMD i więzów (SMD on … constraints yes itp.) była obecna.

template.run    - przykładowy skrypt uruchamiający symulacje na komputerze na jakim 
  zamierzacie liczyć - z podanym wierszem uruchamiającym namd. Pliki wejściowe i wyjściowe będą zdefiniowane jako INPF i OUTF więc tak należy je traktować w wierszu uruchamiającym namda (/home/kasia/NAMD/namd2 +p2 $INPF > $OUTF 2>&)

‘selection’    - wybór atomów zatrzymanych (‘selection constraints’) oraz ciągniętych 
 (‘selection pull’) w symulacji SMD. Jest to zmienna tekstowa, koniecznie w cudzysłowiu ‘’. Konwencja wyboru atomów jest taka jak w MDAnalysis (https://docs.mdanalysis.org/1.1.0/documentation_pages/selections.html) czyli na przykład: ‘protein and segid A B C’  lub  ‘protein and resid 1:55 66:128’. Niestety nie ma ‘chain’ trzeba używać ‘segid’ zamiast tego.
Z podanych zakresów atomów wybrane zostaną atomy CA i do nich zastosujemy opcję constrain lub pull - powstanie plik wejściowy do SMD (SMD_constraints.pdb) zawierający odpowiednie wartości w kolumnach O i B.

Co dostaniemy?
--------------------------------
Program wygeneruje katalog Output zawierający pliki wejściowe, oraz podkatalogi odpowiadające każdemu z kierunków “ciągnięcia”. Każdy taki podkatalog zawiera odpowiednio spreparowany plik wejściowy do namd oraz skrypt bash służący do uruchomienia danej symulacji (na podstawie podanego template.run). Skrypty te uruchamiamy każdy z osobna:

 _. Output/SMD_theta_0_phi_0/run.bash_ (“kropka” uruchomi skrypt tam gdzie jest plik run.bash)
 
lub zbiorczo za pomocą skryptu master.run:

_./master.run_

Uruchomi to symulacje SMD.
