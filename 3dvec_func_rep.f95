module v3d_func_rep
implicit none
contains

  function get_refl_plane(np,u) result(v)
    double precision, dimension(3) :: np, u, v
    
    v(:) = u(:) - 2*get_sc_prod(u,np)/get_sc_prod(np,np)*np(:)
  end function

  function get_refl_plane_n(np,u) result(v)
    double precision, dimension(3) :: np, u, v
    
    v(:) = u(:) - 2*get_sc_prod(u,np)*np(:)
  end function

  function get_rotmat_aa(iaxs,agl) result(rmat)
    double precision :: agl, cosa, sina
    double precision, dimension(3) :: axs, iaxs
    double precision, dimension(3,3) :: rmat
    
    axs = get_normal(iaxs)
  
    cosa = dcos(agl)
    sina = dsin(agl)

    rmat(1:3,1) = (/ axs(1)**2*(1-cosa)+cosa, axs(2)*axs(1)*(1-cosa)+axs(3)*sina, axs(3)*axs(1)*(1-cosa)-axs(2)*sina /)
    rmat(1:3,2) = (/ axs(1)*axs(2)*(1-cosa)-axs(3)*sina, axs(2)**2*(1-cosa)+cosa, axs(3)*axs(2)*(1-cosa)+axs(1)*sina /)
    rmat(1:3,3) = (/ axs(1)*axs(3)*(1-cosa)+axs(2)*sina, axs(2)*axs(3)*(1-cosa)-axs(1)*sina, axs(3)**2*(1-cosa)+cosa /)
  end function
  
  function get_lin_map(mat,u) result(v)
    double precision, dimension(3) :: u, v
    double precision, dimension(3,3) :: mat
    integer :: i, j
    
    do i = 1, 3
      v(i) = 0.d0
      do j = 1, 3
        v(i) = v(i) + mat(i,j)*u(j)
      end do
    end do
  end function
  
  function get_det3x3(mat) result(det)
    double precision, dimension(3,3) :: mat
    double precision :: det

    det = mat(1,1)*mat(2,2)*mat(3,3) + mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(2,1)*mat(3,2)&
    & - mat(1,3)*mat(2,2)*mat(3,1) - mat(1,2)*mat(2,1)*mat(3,3) - mat(1,1)*mat(3,2)*mat(2,3)
  end function

  function get_normal(vec) result(nrmv)
    double precision, dimension(3) :: vec, nrmv
    double precision :: nrm
    nrm = get_norm(vec)
    nrmv(:) = vec(:)/nrm
  end function

  function get_dist_to_line(lp1,lp2,pt) result(dist)
    double precision, dimension(3) :: lp1, lp2, pt, zx, yx, v
    double precision :: dist
    zx = pt - lp1
    yx = get_normal(lp2 - lp1)
    v = zx - get_sc_prod(zx,yx)*yx
    dist = get_norm(v)
  end function
  
  function get_angle(v1,v2) result(agl)
    double precision, dimension(3) :: v1, v2
    double precision :: agl
    
    agl = dacos(get_sc_prod(get_normal(v1),get_normal(v2)))
  end function

  function get_dist(v1,v2) result(dist)
    double precision, dimension(3) :: v1, v2, vec
    double precision :: dist
    vec = v2 - v1
    dist = get_norm(vec)
  end function
  
  function get_norm(vec) result(nrm)
    double precision, dimension(3) :: vec
    double precision :: nrm
    
    nrm = dsqrt(get_sc_prod(vec,vec))
  end function
  
  function get_nrm_sq(vec) result(nrmsq)
    double precision, dimension(3) :: vec
    double precision :: nrmsq
    
    nrmsq = get_sc_prod(vec,vec)
  end function

  function get_sc_prod(v1,v2) result(prod)
    double precision, dimension(3) :: v1, v2
    double precision :: prod
    prod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
  end function

  function get_cr_prod(v1,v2) result(prod)
    double precision, dimension(3) :: v1, v2, prod
    prod(1) = v1(2)*v2(3) - v1(3)*v2(2)
    prod(2) = v1(3)*v2(1) - v1(1)*v2(3)
    prod(3) = v1(1)*v2(2) - v1(2)*v2(1)
  end function
end module
