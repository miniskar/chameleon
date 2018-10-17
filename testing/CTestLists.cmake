#
# Check testing/
#

set(TEST_CMD_shm    testing 4 0 19 7 )
set(TEST_CMD_shmgpu testing 4 1 19 7 )
# set(TEST_CMD_mpi    testing 4 0 19 7 )
# set(TEST_CMD_mpigpu testing 4 1 19 7 )

set( TEST_CATEGORIES shm )
if (CHAMELEON_USE_CUDA AND CUDA_FOUND)
   set( TEST_CATEGORIES ${TEST_CATEGORIES} shmgpu )
endif()

foreach(cat  ${TEST_CATEGORIES})
  foreach(prec ${RP_CHAMELEON_PRECISIONS})

    string(TOUPPER ${prec} PREC)

    if (CHAMELEON_PREC_${PREC})
        add_test(test_${cat}_${prec}lange ./${prec}${TEST_CMD_${cat}} LANGE 117 213 232 )
        add_test(test_${cat}_${prec}gemm  ./${prec}${TEST_CMD_${cat}} GEADD 1.7 -2.3 117 213    215 220    )
        add_test(test_${cat}_${prec}gemm  ./${prec}${TEST_CMD_${cat}} GEMM  1.7 -2.3 117 213 97 215 220 232)
        add_test(test_${cat}_${prec}trsm  ./${prec}${TEST_CMD_${cat}} TRSM      -2.3 117 213    215 220    )
        add_test(test_${cat}_${prec}trmm  ./${prec}${TEST_CMD_${cat}} TRMM      -2.3 117 213    215 220    )
        add_test(test_${cat}_${prec}symm  ./${prec}${TEST_CMD_${cat}} SYMM  1.7 -2.3 117 213    215 220 232)
        add_test(test_${cat}_${prec}syrk  ./${prec}${TEST_CMD_${cat}} SYRK  1.7 -2.3 117 213    215 220    )
        add_test(test_${cat}_${prec}syr2k ./${prec}${TEST_CMD_${cat}} SYR2K 1.7 -2.3 117 213    215 220 232)

        if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
          add_test(test_${cat}_${prec}hemm  ./${prec}${TEST_CMD_${cat}} HEMM  1.7 -2.3 117 213 215 220 232)
          add_test(test_${cat}_${prec}herk  ./${prec}${TEST_CMD_${cat}} HERK  1.7 -2.3 117 213 215 220    )
          add_test(test_${cat}_${prec}her2k ./${prec}${TEST_CMD_${cat}} HER2K 1.7 -2.3 117 213 215 220 232)
        endif()

         add_test(test_${cat}_${prec}posv        ./${prec}${TEST_CMD_${cat}} POSV            117 155 25 143)
         add_test(test_${cat}_${prec}potri       ./${prec}${TEST_CMD_${cat}} POTRI           117 155       )
         add_test(test_${cat}_${prec}gels_qr     ./${prec}${TEST_CMD_${cat}} GELS      0 233 117 255 25 242  )
         add_test(test_${cat}_${prec}gels_hqr    ./${prec}${TEST_CMD_${cat}} GELS      1 233 117 255 25 242 3)
         add_test(test_${cat}_${prec}gels_lq     ./${prec}${TEST_CMD_${cat}} GELS      0 117 233 155 25 242  )
         add_test(test_${cat}_${prec}gels_hlq    ./${prec}${TEST_CMD_${cat}} GELS      1 117 233 155 25 242 3)
         add_test(test_${cat}_${prec}gesv_incpiv ./${prec}${TEST_CMD_${cat}} GESV_INCPIV     117 155 25 242)

         add_test(test_${cat}_${prec}gels_hqr_greedy    ./${prec}${TEST_CMD_${cat}} GELS_HQR 233 117 255 25 242 2 -1 1 -1 0)
         add_test(test_${cat}_${prec}gels_hqr_fibonacci ./${prec}${TEST_CMD_${cat}} GELS_HQR 233 117 255 25 242 2 -1 2 -1 0)
         add_test(test_${cat}_${prec}gels_hqr_binary    ./${prec}${TEST_CMD_${cat}} GELS_HQR 233 117 255 25 242 2 -1 3 -1 0)
         add_test(test_${cat}_${prec}gels_hlq_greedy    ./${prec}${TEST_CMD_${cat}} GELS_HQR 117 233 255 25 242 2 -1 1 -1 0)
         add_test(test_${cat}_${prec}gels_hlq_fibonacci ./${prec}${TEST_CMD_${cat}} GELS_HQR 117 233 255 25 242 2 -1 2 -1 0)
         add_test(test_${cat}_${prec}gels_hlq_binary    ./${prec}${TEST_CMD_${cat}} GELS_HQR 117 233 255 25 242 2 -1 3 -1 0)

         add_test(test_${cat}_${prec}gels_qr_systolic   ./${prec}${TEST_CMD_${cat}} GELS_SYSTOLIC 233 117 255 25 242 3 2)
         add_test(test_${cat}_${prec}gels_lq_systolic   ./${prec}${TEST_CMD_${cat}} GELS_SYSTOLIC 117 233 255 25 242 3 2)
     endif()

  endforeach()
endforeach()

#if (CHAMELEON_USE_MPI AND MPI_C_FOUND)
#    set( TEST_CATEGORIES ${TEST_CATEGORIES} mpi )
#    if (CHAMELEON_USE_CUDA AND CUDA_FOUND)
#        set( TEST_CATEGORIES ${TEST_CATEGORIES} mpigpu )
#    endif()
#    foreach(prec ${RP_CHAMELEON_PRECISIONS})
#        add_test(test_mpi_${prec}lange mpirun -np 4 ./${prec}testing 1 0 LANGE 600 500 600 --p=2)
#    endforeach()
#endif()
