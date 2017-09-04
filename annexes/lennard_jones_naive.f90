do io=1,no

    rotx=[grid%rotxx(io),grid%rotxy(io),grid%rotxz(io)]
    roty=[grid%rotyx(io),grid%rotyy(io),grid%rotyz(io)]
    rotz=[grid%rotzx(io),grid%rotzy(io),grid%rotzz(io)]

    do ix=1,nx
      do iy=1,ny
        do iz=1,nz

          ! pour chaque solvant
          ! on calcule la position de chaque site du solvant
          do s=1,ns
            do v=1,size(solvent(s)%site)
              do u=1,size(solute%site)
                epsuv(u,v)=sqrt(solute%site(u)%eps * solvent(s)%site(v)%eps)
                siguv(u,v)=(solute%site(u)%sig + solvent(s)%site(v)%sig)/2._dp
              end do
            end do
            vloc=0._dp
            do ss=1,size(solvent(s)%site)
              if( solvent(s)%site(ss)%eps<=epsdp) cycle
              xss=x(ix)+dot_product(rotx,solvent(s)%site(ss)%r)
              yss=y(iy)+dot_product(roty,solvent(s)%site(ss)%r)
              zss=z(iz)+dot_product(rotz,solvent(s)%site(ss)%r)
              ! pour chaque site de solute
              ! on calcule la distance entre le site de solute et le site de solvant
              ! on calcule vlj(epsij,sigij,rsq)
              ! et on ajoute la contribution a v
              do u=1,size(solute%site)
                if( solute%site(u)%eps <= epsdp .or. vloc > 1.e5 ) cycle ! if the solute site does not wear a Lennard-Jones contribution
                dx =abs(xss-solute%site(u)%r(1)); do while(dx>lx/2._dp); dx=abs(dx-lx); end do
                dy =abs(yss-solute%site(u)%r(2)); do while(dy>ly/2._dp); dy=abs(dy-ly); end do
                dz =abs(zss-solute%site(u)%r(3)); do while(dz>lz/2._dp); dz=abs(dz-lz); end do
                vloc = vloc + vlj( epsuv(u,ss), siguv(u,ss), dx**2+dy**2+dz**2)
              end do
            end do

            solvent(s)%vext(io,ix,iy,iz) = solvent(s)%vext(io,ix,iy,iz) + vloc

          end do

        end do
      end do
    end do

  end do