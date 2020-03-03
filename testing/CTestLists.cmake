#
# Check testing/
#
set(NP 2) # Amount of MPI processes
set(THREADS 2) # Amount of threads
set(N_GPUS 0) # Amount of graphic cards
set(TEST_CATEGORIES shm)
if (CHAMELEON_USE_MPI AND MPI_C_FOUND)
  set( TEST_CATEGORIES ${TEST_CATEGORIES} mpi )
endif()
if (CHAMELEON_USE_CUDA AND CUDA_FOUND)
  set(N_GPUS 0 1)
endif()

foreach(prec ${RP_CHAMELEON_PRECISIONS})
  set (CMD ./chameleon_${prec}testing)

  #
  # Create the list of test based on precision and runtime
  #
  set( TESTS lacpy lange lantr lansy )
  if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
    set( TESTS ${TESTS} lanhe )
  endif()
  set( TESTS ${TESTS}
    geadd tradd lascal
    gemm symm syrk syr2k trmm trsm )
  if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
    set( TESTS ${TESTS}
      hemm herk her2k )
  endif()
  set( TESTS ${TESTS}
    potrf potrs posv trtri lauum )
  if ( NOT CHAMELEON_SCHED_PARSEC )
    set( TESTS ${TESTS} potri )
  endif()
  if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
    set( TESTS ${TESTS}
      sytrf sytrs sysv )
  endif()
  set( TESTS ${TESTS}
    getrf     getrs gesv
    geqrf     gelqf
    geqrf_hqr gelqf_hqr)
  if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
    set( TESTS ${TESTS}
      ungqr     unglq     unmqr     unmlq
      ungqr_hqr unglq_hqr unmqr_hqr unmlq_hqr)
  else()
    set( TESTS ${TESTS}
      orgqr     orglq     ormqr     ormlq
      orgqr_hqr orglq_hqr ormqr_hqr ormlq_hqr)
  endif()
  set( TESTS ${TESTS}
    #geqrs     gelqs
    #geqrs_hqr gelqs_hqr
    gels
    gels_hqr )

  foreach(cat ${TEST_CATEGORIES})
    foreach(gpus ${N_GPUS})

      if (${gpus} EQUAL 1)
        set(cat ${cat}_gpu)
      endif()

      if (${cat} STREQUAL "mpi")
        set (PREFIX mpiexec --bind-to none -n ${NP})
      else()
        set (PREFIX "")
      endif()

      foreach(_test ${TESTS})
        add_test(test_${cat}_${prec}${_test} ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/${_test}.in )
      endforeach()
    endforeach()
  endforeach()
endforeach()
