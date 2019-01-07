%
% RANDOM GAMMA NUMBER GENERATOR 
%
% IT GENERATES RANDOM NUMBERS FROM A GAMMA DISTRIBUTION 
% WITH PDF  Y^(A-1) EXP(-A*Y) UP TO A NORMALIZING CONSTANT
% WITH PARAMETERS A, B AND MEAN E(X|A,B) = A/B
%
% =========================================================================
%
% CALL :
%
% Y = rngGamma( A ) , A > 0
% Y = rngGamma( A , B ), A > 0, B > 0
%
% =========================================================================
%
% REFERENCES :
%
% Ahrens, J. H. and Dieter, U. (1982). 
% Generating gamma variates by a modified rejection technique. 
% Communications of the ACM, 25, 47–54. Algorithm GD, p. 53.
%
% Ripley, BD (1987). 
% Stochastic simulation. 
% New York: Wiley. 
%
% Marsaglia, G., and W. W. Tsang. (2000).
% "A Simple Method for Generating Gamma Variables." 
% ACM Transactions on Mathematical Software. Vol. 26, p. 363–372. 

function y = rngGamma(a,b) 

    if nargin > 2
        error(message('rngGamma:TooManyInputs')) ;
    end

    if a > 1
       y = rngGamma1(a) ; 
    elseif a == 1
        y = -log(rand()) ;
    elseif 0<a && a < 1
       y = rngGamma2(a) ; 
    else
        error(message('rngGamma: Invalid value a <= 0')) ;
    end
    
    if nargin == 2
        if b<=0 
            error(message('rngGamma: Invalid value b <= 0')) ;
        end
        y = y/b ;
    end
    
end

% =========================================================================

function y = rngGamma1(s)

    d = s - 0.333333333333333 ;
    c = 1.0/sqrt(9.0*d) ;

    while true

      while true
        x = randn() ;
        v = (1.0 + c*x)^3 ;
        if v > 0.0 ; break ; end
      end

      u = randn() ;
      if (u < 1.0 - 0.0331*x^4)
        y = d*v ;
        break
      elseif log(u) < 0.5*x^2 + d*(1.0-v+log(v))
        y = d*v ;
        break
      end

    end

end

% =========================================================================

function y = rngGamma2(s)

    e = 2.718281828459046 ;

    b = (s+e)/e ;
    c1 = 1.0/s ;
    
    while true
        bu = b*rand() ;
		if bu <= 1.0
 		    y = exp(max(-40,c1*log(bu))) ;
		    if rand() < exp(-y) 
                return
            end
		else
		    y = -log((b-bu)/s) ;
		    if rand() < y^(s-1)
                return
            end
        end
    end
    
end

% =========================================================================

% TO BE DONE :
%
% UPDATE rngGamma1 using newer algorithm
% TIDY UP rngGamma2

% =========================================================================

