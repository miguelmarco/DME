all: encrypt decrypt encrypt_with_skey skey_to_pkey test kat-nist kat

hash.o: hash.c hash.h
	gcc -c hash.c -Wall -Wextra -std=c99 -pedantic

bigint.o: bigint.c bigint.h
	gcc -c bigint.c -Wall -Wextra -std=c99 -pedantic

dme.o: dme.c dme.h setup.h hash.h
	gcc -c dme.c -Wall -Wextra -std=c99 -pedantic

utils.o: utils.h utils.c
	gcc -c utils.c -Wall -Wextra -std=c99 -pedantic

encrypt: encrypt.c dme.o utils.o bigint.o hash.o
	gcc -o encrypt encrypt.c dme.o utils.o bigint.o hash.o -Wall -Wextra -std=c99 -pedantic

decrypt: decrypt.c dme.o utils.o bigint.o hash.o
	gcc -o decrypt decrypt.c dme.o utils.o bigint.o hash.o -Wall -Wextra -std=c99 -pedantic

encrypt_with_skey: encrypt_with_skey.c dme.o utils.o bigint.o hash.o
	gcc -o encrypt_with_skey encrypt_with_skey.c dme.o utils.o bigint.o hash.o -Wall -Wextra -std=c99 -pedantic

skey_to_pkey: skey_to_pkey.c dme.o utils.o bigint.o hash.o
	gcc -o skey_to_pkey skey_to_pkey.c dme.o utils.o bigint.o hash.o -Wall -Wextra -std=c99 -pedantic

test: test.c dme.o utils.o bigint.o hash.o
	gcc -o test test.c dme.o utils.o bigint.o hash.o -Wall -Wextra -std=c99 -pedantic -lm

kat: test
	./test 500 > kat.txt

find_min_poly: find_min_poly.c
	gcc -o find_min_poly find_min_poly.c -Wall -Wextra -std=c99 -pedantic
	
find_best_exps: find_best_exps.c bigint.c bigint.h
	gcc -o find_best_exps find_best_exps.c bigint.c -Wall -Wextra -std=c99 -pedantic

setup1.h: find_min_poly find_best_exps
	./find_min_poly > setup1.h
	./find_best_exps >> setup1.h

find_m1_m2: find_m1_m2.c setup1.h
	gcc -o find_m1_m2 find_m1_m2.c -Wall -Wextra -std=c99 -pedantic

setup.h: setup1.h find_m1_m2
	cat setup1.h > setup.h
	./find_m1_m2 >> setup.h

dme-nist.o: dme.c dme.h rng.h api.h setup.h
	gcc -o dme-nist.o -c dme.c -DNIST -Wall -Wextra -std=c99 -pedantic

rng.o: rng.h
	gcc -c rng.c

PQCgenKAT_kem.o: PQCgenKAT_kem.c rng.h api.h
	gcc -c PQCgenKAT_kem.c -Wall -Wextra -std=c99 -pedantic

PQCgenKAT_kem: PQCgenKAT_kem.o dme-nist.o rng.o hash.o
	gcc -o PQCgenKAT_kem PQCgenKAT_kem.o dme-nist.o rng.o hash.o -Wall -Wextra -lcrypto -pedantic

kat-nist: PQCgenKAT_kem
	./PQCgenKAT_kem

clean:
	rm -f setup.h setup1.h
	rm -f *.o
	rm -f encrypt decrypt encrypt_with_skey skey_to_pkey test
	rm -f kat.txt stats.txt
	rm -f find_min_poly find_best_exps find_m1_m2
	rm -f PQCgenKAT_kem
	rm -f *.req *.rsp

.PHONY: clean all kat kat-nist
