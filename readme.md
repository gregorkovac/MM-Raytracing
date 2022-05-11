## Todo in razne opombe

**Funkcija `raytracing`**:
  - [x] Iskanje presečišča: Gregor
  - [ ] Odbojni kot: Lana
  - [ ] Barvanje: Jovana

**To fix:** 
- [ ] Zaenkrat so parcialni odvodi ravnine argumenti funkcije. Treba je pogruntat, če se jih
da izračunat v Octave.
- [ ] Računanje ray hita je mogoče mal ugly

**Testiranje funkcije:** \
Napišeš funkcijo ravnine `f`, napišeš odvode po vsaki spremenljivki (`dfdx`, `dfdy` in `dfdz`), ustvariš začetno točko raytracinga `T0`, ustvariš smerni vektor `v`
in pokličeš funkcijo.
```
f = @(x, y, z) x+y+z;
dfdx = @(x, y, z) 1
dfdy = @(x, y, z) 1
dfdz = @(x, y, z) 1
T0 = [5; 5; 5];
v = [-1; -1; -1];
raytracing(f, dfdx, dfdy, dfdz, T0, v)
```