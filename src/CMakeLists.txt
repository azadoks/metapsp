###########################################################
# METAPSP library
###########################################################
set(src_libmetapsp
    aeo.f90
    check_data.f90
    const_basis.f90
    derivs.f90
    dpnint.f90
    eresid.f90
    exc_libxc_stub.f90
    exc_off_pbe.f
    exchdl.f90
    excpzca.f90
    excr2scan01_m.f90
    excwig.f90
    fcrossover.f90
    fphsft_m.f90
    fpovlp.f90
    gauss_smooth.f90
    globalmath_loc.F90
    gnu_script_m.f90
    linout.f90
    lschfb.f90
    lschfb_m.f90
    lschfs_m.f90
    lschkb_m.f90
    lschps_m.f90
    lschpsbar_m.f90
    lschvkbb_m.f90
    lschvkbbe_m.f90
    lschvkbs_m.f90
    optimize.f90
    polysmooth.f90
    pspot_m.f90
    qroots.f90
    r2scanmod.F90
    rtloc_m.f90
    run_diag_m.f90
    run_ghosts_m.f90
    run_optimize_m.f90
    run_phsft_m.f90
    run_plot_m.f90
    run_vkb_m.f90
    sbf8.f90
    sbf_basis.f90
    sbf_basis_con.f90
    sbf_rc_der.f90
    sratom.f90
    sratom_m.f90
    tfapot.f90
    upfout_m.f90
    vkboutwf_m.f90
    vkbphsft_m.f90
    vout.f90
    vout_m.f90
    vpinteg.f90
    vploc.f90
    wellstate_m.f90
    wf_rc_der.f90
)
if (WITH_PSML)
    list(APPEND src_libmetapsp m_psmlout_m.f90 m_getopts.f90 m_uuid.f90)
endif()
add_library(metapsp STATIC ${src_libmetapsp})
add_dependencies(metapsp)
target_link_libraries(metapsp ${COMMON_LIBRARIES})

###########################################################
# oncvpspm.x (with PSML support)
# oncvpspmx.x (without PSML support)
###########################################################
if (WITH_PSML)
    add_executable(oncvpspm oncvpsp_m.f90)
    set_target_properties(oncvpspm
        PROPERTIES
            OUTPUT_NAME oncvpspm.x
    )
    add_dependencies(oncvpspm metapsp)
    target_link_libraries(oncvpspm metapsp)
else()
    add_executable(oncvpspmx oncvpsp_mx.f90)
    set_target_properties(oncvpspmx
        PROPERTIES
            OUTPUT_NAME oncvpspmx.x
    )
    add_dependencies(oncvpspmx metapsp)
    target_link_libraries(oncvpspmx metapsp)
endif()
