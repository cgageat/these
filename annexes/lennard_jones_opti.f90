       do s=1,ns
          do ss=1,solvent(s)%nsite
            do io=1,no
              xmod (io,ss,s) = DOT_PRODUCT( [grid%rotxx(io), grid%rotxy(io), grid%rotxz(io)] , solvent(s)%site(ss)%r )
              ymod (io,ss,s) = DOT_PRODUCT( [grid%rotyx(io), grid%rotyy(io), grid%rotyz(io)] , solvent(s)%site(ss)%r )
              zmod (io,ss,s) = DOT_PRODUCT( [grid%rotzx(io), grid%rotzy(io), grid%rotzz(io)] , solvent(s)%site(ss)%r )
            end do
          end do
        end do





        do u=1,size(solute%site)
          if( solute%site(u)%eps <= epsdp ) cycle ! if the solute site does not wear a Lennard-Jones contribution

          ! prepare table to loop over
          ! These tables allow to loop only on a cube with l=cutoff
          minx = floor(mod((solute%site(u)%r(1)-cutoff+lx), lx)/grid%dx)
          maxx = floor(mod((solute%site(u)%r(1)+cutoff   ), lx)/grid%dx)
          miny = floor(mod((solute%site(u)%r(2)-cutoff+ly), ly)/grid%dy)
          maxy = floor(mod((solute%site(u)%r(2)+cutoff   ), ly)/grid%dy)
          minz = floor(mod((solute%site(u)%r(3)-cutoff+lz), lz)/grid%dz)
          maxz = floor(mod((solute%site(u)%r(3)+cutoff   ), lz)/grid%dz)

          if (minx<maxx) then
            xtab = (/ (I, I = minx, maxx) /)
          else ! take into account PBC
            xtab = [(/ (I, I = 0, maxx) /), (/ (I, I = minx, nx) /)]
          end if

          if (miny<maxy) then
            ytab = (/ (I, I = miny, maxy) /)
          else ! take into account PBC
            ytab = [(/ (I, I = 0, maxy) /), (/ (I, I = miny, ny) /)]
          end if

          if (minz<maxz) then
            ztab = (/ (I, I = minz, maxz) /)
          else ! take into account PBC
            ztab = [(/ (I, I = 0, maxz) /), (/ (I, I = minz, nz) /)]
          end if
          !end of prepare table to loop over


          do s=1,ns
            do ss=1,size(solvent(s)%site)
              if( solvent(s)%site(ss)%eps<=epsdp ) cycle
              epsuv=sqrt(solute%site(u)%eps * solvent(s)%site(ss)%eps)
              siguv6=(  (solute%site(u)%sig + solvent(s)%site(ss)%sig)/2._dp)**6

              do indextabz=1,ztabsize ! indextabz is the index in ztabsize
                iz = ztab(indextabz)  ! iz is the index of the point in the grid
                zgrid=z(iz)           ! zgrid is the z position of the point

                do indextaby=1,ytabsize
                  iy = ytab(indextaby)
                  ygrid=y(iy)

                  do indextabx=1,xtabsize
                    ix = xtab(indextabx)
                    xgrid=x(ix)

                    do io=1,no 
                      if( solvent(s)%vext(io,ix,iy,iz) > 1.e5 ) cycle

                      xss=xgrid+xmod(io,ss,s)
                      yss=ygrid+ymod(io,ss,s)
                      zss=zgrid+zmod(io,ss,s)

                      dx =abs(xss-solute%site(u)%r(1)); do while(dx>lx/2._dp); dx=abs(dx-lx); end do
                      dy =abs(yss-solute%site(u)%r(2)); do while(dy>ly/2._dp); dy=abs(dy-ly); end do
                      dz =abs(zss-solute%site(u)%r(3)); do while(dz>lz/2._dp); dz=abs(dz-lz); end do

                      solvent(s)%vext(io,ix,iy,iz) = solvent(s)%vext(io,ix,iy,iz) + vlj( epsuv, siguv6, dx**2+dy**2+dz**2, cutoffsq)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do