module chion
    ! The chion facade: the only module a host program needs to `use`.
    !
    ! Mirrors yelmo/src/yelmo.f90 -- a bare re-export wrapper with no
    ! declarations and no code of its own. `use` without `only` re-exports the
    ! public entities of each module, so `use chion` brings in the public API,
    ! the derived types, the precision kinds and the model state types
    ! together.
    !
    ! Typical host usage:
    !
    !     use chion
    !
    !     type(chion_class) :: chn
    !     real(wp) :: smb(ncol)
    !
    !     call chion_init(chn,"par/chion_column.nml",ncol)
    !     call chion_init_state(chn)
    !     do n = 1, nstep
    !         chn%forc%air_temperature(:) = ...
    !         chn%forc%snowfall_rate(:)   = ...
    !         call chion_update(chn,dt_days)
    !         call chion_get_smb(chn,smb)          ! [kg m-2 s-1], + = ice gains
    !     end do
    !     call chion_end(chn)
    !
    ! The model state types are re-exported so that a host or a diagnostic can
    ! reach into chn%bsi / chn%pdd / chn%itm for anything the model-agnostic
    ! API does not expose (layer profiles, runoff, albedo, ...). Only the
    ! state of the SELECTED model is allocated.

    use chion_defs
    use chion_model
    use chion_api
    use chion_io

    use snow_bessi
    use snow_pdd
    use snow_itm

    use snow_diagnostics

    use chion_forcing_monthly

end module chion
