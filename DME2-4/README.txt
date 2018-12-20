/

    README_DME-(4-2-48)  : this file


/Reference_Implementation/

    api.h           : the required api
    dme.h and dme.c : the implementation of the DME primitives
    rng.h and rng.c : random number generator (provided by NIST)
    PQCgenKAT_kem.c : the known answer test generator (provided by NIST)
    makefile        : for automated building (setup.h, dme.o, KAT

    bigint.h and bigint.c : fixed precision 384-bits integer arithmetic
    find_min_poly.c       : to find random irreducible polynomials of degree 2
                            and 3 in Fq[X]
    find_best_exps.c      : to find the best exponential maps and related
                            constants
    find_m1_m2.c          : to precompute the matrices M1, M2, M1^(-1), M2^(-1)
                            used to convert a secret key into a public key
     hash.h and  hash.c   :  the implementation of the hash-256
    
     decrypt.c              :  use: ./decrypt cipher_text secret_key   [padding]
     encrypt.c             :   use: ./encrypt plain_text public_key   [padding]              
     encrypt_with_skey.c   :  use: ./encrypt_with_skey plain_text secret_key [padding]
               
     skey_to_pkey.c        :  use: ./skey_to_pkey secret_key
    utils.h and utils.c    
   test.c     : Time using for DME-(4,2,48)  using our own APIs.
 
