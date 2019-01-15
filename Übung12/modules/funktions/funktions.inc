function sinus(f) result(r)
    type(funktion) :: f,r,s
    procedure(func),pointer :: p_sin => p_sinus
    s = create_funktion(p_sin)
    r = base_chain(s,f)
end function

real function p_sinus(x) 
     real ::x 
     p_sinus = sin(x)
end function

function cosinus(f) result(r)
    type(funktion) :: f,r,s
    procedure(func),pointer :: p_cosin => p_cosinus
    s = create_funktion(p_cosin)
    r = base_chain(s,f)
end function

real function p_cosinus(x) 
     real ::x 
     p_cosinus = cos(x)
end function


function euler(f) result(r)
    type(funktion) :: f,r,s
    procedure(func),pointer :: p_exp => p_euler
    s = create_funktion(p_exp)
    r = base_chain(s,f)
end function

real function p_euler(x) 
     real ::x 
     p_euler = exp(x)
end function

function heaviside(f) result(r)
    type(funktion) :: f,r,s
    procedure(func),pointer :: p_hs => p_heaviside
    s = create_funktion(p_hs)
    r = base_chain(s,f)
end function

real function p_heaviside(x) 
     real ::x 
     if (x>=0) then
       p_heaviside=1.
     else
        p_heaviside=0.
     end if
end function



