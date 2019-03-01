#
# Check timing/
#

set(TEST_CMD_shm    --n_range=17:407:39    --nb=32  --ib=7                )
set(TEST_CMD_shmgpu --n_range=170:4070:390 --nb=320 --ib=48 --gpus=1      )
set(TEST_CMD_mpi    --n_range=17:407:39    --nb=32  --ib=7  --P=2         )
set(TEST_CMD_mpigpu --n_range=170:4070:390 --nb=320 --ib=48 --P=2 --gpus=1)

set(MPI_CMD_shm    )
set(MPI_CMD_shmgpu )
set(MPI_CMD_mpi    mpirun -np 4)
set(MPI_CMD_mpigpu mpirun -np 4)

if (NOT CHAMELEON_SIMULATION)
  set( TEST_CATEGORIES shm )
  if (CHAMELEON_USE_CUDA AND CUDA_FOUND)
    set( TEST_CATEGORIES ${TEST_CATEGORIES} shmgpu )
  endif()
else()
  set( TEST_CATEGORIES simushm )
  if (CHAMELEON_USE_CUDA AND CUDA_FOUND)
    set( TEST_CATEGORIES ${TEST_CATEGORIES} simugpu )
  endif()
endif()

set(TESTLIST
  gels
  gemm
  getrf_incpiv
  getrf_nopiv
  geqrf
  gelqf
  posv
  potrf
  potri
  )

set(CHAMELEON_PRECISIONS_ZC "c;z")
set(TESTLIST_ZC
  sytrf
  )

if (NOT CHAMELEON_SIMULATION)

  foreach(cat ${TEST_CATEGORIES})
    foreach(prec ${RP_CHAMELEON_PRECISIONS})
      string(TOUPPER ${prec} PREC)
      if (CHAMELEON_PREC_${PREC})
        foreach(test ${TESTLIST})
          add_test(time_${cat}_${prec}${test} ${MPI_CMD_${cat}} ./time_${prec}${test}_tile ${TEST_CMD_${cat}} --check)
        endforeach()
      endif()
    endforeach()
    foreach(prec ${CHAMELEON_PRECISIONS_ZC})
      string(TOUPPER ${prec} PREC)
      if (CHAMELEON_PREC_${PREC})
        foreach(test ${TESTLIST_ZC})
          add_test(time_${cat}_${prec}${test} ${MPI_CMD_${cat}} ./time_${prec}${test}_tile ${TEST_CMD_${cat}} --check)
        endforeach()
      endif()
    endforeach()
  endforeach()

  if ( CHAMELEON_SCHED_STARPU )
    foreach(cat ${TEST_CATEGORIES})
      foreach(prec ${RP_CHAMELEON_PRECISIONS})
        string(TOUPPER ${prec} PREC)
        if (CHAMELEON_PREC_${PREC})
          foreach(test ${TESTLIST})
            add_test(time_ooc_${cat}_${prec}${test} ${MPI_CMD_${cat}} ./time_${prec}${test}_tile ${TEST_CMD_${cat}} --ooc )#--check)
            set_tests_properties(time_ooc_${cat}_${prec}${test} PROPERTIES
              ENVIRONMENT "STARPU_DISK_SWAP=/tmp;STARPU_LIMIT_CPU_MEM=1"
              )
          endforeach()
        endif()
      endforeach()
      foreach(prec ${CHAMELEON_PRECISIONS_ZC})
        string(TOUPPER ${prec} PREC)
        if (CHAMELEON_PREC_${PREC})
          foreach(test ${TESTLIST_ZC})
            add_test(time_ooc_${cat}_${prec}${test} ${MPI_CMD_${cat}} ./time_${prec}${test}_tile ${TEST_CMD_${cat}} --ooc )#--check)
            set_tests_properties(time_ooc_${cat}_${prec}${test} PROPERTIES
              ENVIRONMENT "STARPU_DISK_SWAP=/tmp;STARPU_LIMIT_CPU_MEM=1"
              )
          endforeach()
        endif()
      endforeach()
    endforeach()
  endif()

  if (CHAMELEON_USE_MPI AND MPI_C_FOUND)
    set( TEST_CATEGORIES mpi )
    #set( TEST_CATEGORIES ${TEST_CATEGORIES} mpi )
    #if (CHAMELEON_USE_CUDA AND CUDA_FOUND)
    #    set( TEST_CATEGORIES ${TEST_CATEGORIES} mpigpu )
    #endif()
    set(TESTLIST_MPI
      potrf
      )
    foreach(cat ${TEST_CATEGORIES})
      foreach(prec ${RP_CHAMELEON_PRECISIONS})
        string(TOUPPER ${prec} PREC)

        if (CHAMELEON_PREC_${PREC})
          foreach(test ${TESTLIST_MPI})
            add_test(time_${cat}_${prec}${test} ${MPI_CMD_${cat}} ./time_${prec}${test}_tile ${TEST_CMD_${cat}} --check)
          endforeach()
        endif()
      endforeach()
    endforeach()
  endif()

else (NOT CHAMELEON_SIMULATION)

  set(TEST_CMD_simushm -N 9600:9600:1 -b 960)
  set(TEST_CMD_simugpu -N 9600:9600:1 -b 960 -g 1)
  set(RP_CHAMELEON_PRECISIONS_SIMU "s;d")
  foreach(cat ${TEST_CATEGORIES})
    foreach(prec ${RP_CHAMELEON_PRECISIONS_SIMU})
      string(TOUPPER ${prec} PREC)
      if (CHAMELEON_PREC_${PREC})
        add_test(time_${cat}_${prec}potrf ${MPI_CMD_${cat}} ./time_${prec}potrf_tile ${TEST_CMD_${cat}})
        set_tests_properties(time_${cat}_${prec}potrf PROPERTIES
                             ENVIRONMENT "STARPU_HOME=${CMAKE_SOURCE_DIR}/simucore/perfmodels;STARPU_HOSTNAME=sirocco"
                             )
      endif()
    endforeach()
  endforeach()

endif (NOT CHAMELEON_SIMULATION)
