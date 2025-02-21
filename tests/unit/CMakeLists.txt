set(UNIT_TESTS_NAMES blocs_static dummy indexdw operators operators_static utils_random_matrix
			 simple_blas_add
			 simple_blas_dnrm
			 simple_blas_mul_dgemm
		     simple_blas_mul_dsymm
		     simple_blas_mul_zgemm
		     simple_blas_mul_zhemm
		     simple_lapack_svd_dgesdd
		     simple_lapack_svd_dgesvd
		     simple_lapack_svd_zgesvd
		     # simple_lapack_svd_zgesdd
		     # simple_lapack_eig_zheev
			 routine_svd_routine
			 routine_eig_lanczos
			 routine_mul_routine
			 routine_indices
			 routine_interface
			 qnum
 )

foreach(UNIT_TEST_NAME ${UNIT_TESTS_NAMES})
  add_executable(test_${UNIT_TEST_NAME} test_${UNIT_TEST_NAME}.cpp)
  target_compile_definitions(
    test_${UNIT_TEST_NAME} PUBLIC BOOST_TEST_MODULE=${UNIT_TEST_NAME}
                                  BOOST_TEST_DYN_LINK=TRUE)
  target_link_libraries(
    test_${UNIT_TEST_NAME} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${TBB_TBB_LIBRARY_RELEASE} lapack)
  add_test(unit_${UNIT_TEST_NAME} test_${UNIT_TEST_NAME})
endforeach()
