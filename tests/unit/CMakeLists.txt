set(UNIT_TESTS_NAMES blocs_static dummy indexdw operators operators_static utils_random_matrix
			 blas_add
			 blas_dnrm
			 blas_mul_dgemm
		     blas_mul_dsymm
		     blas_mul_zgemm
		     blas_mul_zhemm
		     lapack_svd_dgesdd
		     lapack_svd_dgesvd
		     lapack_svd_zgesvd
		     lapack_svd_zgesdd
		     # lapack_eig_zheev
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
