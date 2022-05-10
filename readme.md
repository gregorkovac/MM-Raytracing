## Todo in razne opombe

**Funkcija `raytracing`**:
  - [ ] Iskanje presečišča: Gregor
  - [ ] Odbojni kot: Lana
  - [ ] Barvanje: Jovana

Funkcija zaenkrat poišče ray hit, ampak je treba vse skupaj naredit malo bolj natančno z Newtonovo iteracijo.

**Testiranje funkcije:** \
Napišeš funkcijo ravnine `f`, ustvariš začetno točko raytracinga `T0`, ustvariš smerni vektor `v`
in pokličeš funkcijo.
```
f = @(x, y, z) x+y+z;
T0 = [2; 2; 2];
v = [-1; -1; -1];
raytracing(f, T0, v)
```