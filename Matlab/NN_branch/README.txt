2D -implementaatio NN-algoritmiikasta.
Funktioita muokattu yhteensopiviksi ja muutenkin asioita muuteltu.

T�m� olisi hyvin ik�v� sulattaa varsinaiseen versioon, esim.
ett� naapurienetsint�algon voisi valita binaarimuuttujalla.
Mm. random size(blaa,1):ta on muuteltu size(blaa, 2) yms.
Siksi erillinen branch.

Tuon NN-laskennan nopeus on aika mutkikas juttu. 

h_grid = etsint�laatikoiden leveys (testcase.m)
range = etsint�s�de (neighbors_splinesNN.m)

Molemmat on atm vakioita, eiv�t skaalaa kuten teid�n hieno h.
Jos jompaa kumpaa skaalaisi, niin rangea hiukkaskohtaisesti.

Neli�r�j�hde -case ei nopeudu NN:ll�, koska hiukkaset vaikuttavat
toisiinsa liian vahvasti: pit�� joka tapauksessa loopata kaikkien
yli, jotta saa hyvi� tuloksia. NN antaa dramaattisia nopeutuksia,
kun mallinnetaan esim. pitki� tankoja.

Sainhan min� t�m�n lopulta toimimaan.
Saa toki tarkistaa toimiiko teid�nkin mielest�.

