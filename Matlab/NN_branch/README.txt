2D -implementaatio NN-algoritmiikasta.
Funktioita muokattu yhteensopiviksi ja muutenkin asioita muuteltu.

Tämä olisi hyvin ikävä sulattaa varsinaiseen versioon, esim.
että naapurienetsintäalgon voisi valita binaarimuuttujalla.
Mm. random size(blaa,1):ta on muuteltu size(blaa, 2) yms.
Siksi erillinen branch.

Tuon NN-laskennan nopeus on aika mutkikas juttu. 

h_grid = etsintälaatikoiden leveys (testcase.m)
range = etsintäsäde (neighbors_splinesNN.m)

Molemmat on atm vakioita, eivät skaalaa kuten teidän hieno h.
Jos jompaa kumpaa skaalaisi, niin rangea hiukkaskohtaisesti.

Neliöräjähde -case ei nopeudu NN:llä, koska hiukkaset vaikuttavat
toisiinsa liian vahvasti: pitää joka tapauksessa loopata kaikkien
yli, jotta saa hyviä tuloksia. NN antaa dramaattisia nopeutuksia,
kun mallinnetaan esim. pitkiä tankoja.

Sainhan minä tämän lopulta toimimaan.
Saa toki tarkistaa toimiiko teidänkin mielestä.

