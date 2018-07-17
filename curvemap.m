# curvemap.m - Graphically map a low-prime EC curve:
#
# Assumptions: using Weirstrauss (?) form y^2 = x^3 + a*x + b
#
# SOME INTERESTING CURVES:
#
#  (p, a, b, Gidx) = (23,  20,  4, 1);  ==>  (N, n, h) = (30, 30);
#                    (23,  20,  4, 1);  ==>  (N, n, h) = (30, 30);
#                    (23,  20,  4, 1);  ==>  (N, n, h) = (30, 30);
#                    (43,  40,  6, 1);  ==>  (N, n, h) = (50, 25);
#                    (43,  40, 10, 1);  ==>  (N, n, h) = (37, 37);
#                    (43,  40, 12, 1);  ==>  (N, n, h) = (47, 47);
#                   (255, 124, 12, 1);  ==>  (N, n, h) = (121, 2, 60.5);  non-prime field; many roots on some x's
#
# QUESTIONS:
#
#   * Do I include the point at infinity when counting curve points for n?
#
#   * How to represent the point at infinity?
#
# LINKS:
#
#  "Why do public keys need to be validated?"
#   https://crypto.stackexchange.com/questions/3820/why-do-public-keys-need-to-be-validated
#
function curvemap(p, a, b,
                  GenIndex = 1   # Use >1 to pick a later point for generator
                 )
  OrderGuess = 2*p;

  Warnings = "";
  if (length(factor(p))>1)
    if (length(unique(factor(p)))!=1)
      Warnings = sprintf("%s(Notice: Field order p = %g is not a prime power. Field theorems do not hold.)\n", Warnings, p);
    else
      Warnings = sprintf("%s(Notice: Field order p = %g is not prime.)\n", Warnings, p);
    end
  end
  printf(Warnings);

  global param;
  global X;
  global Xp2;
  global Pinf;
  param = [p, a, b];
  X = [0:(p-1)];
  Xp2 = mod(X.^2,p);    # Precompute x^2

  N = 1;          # Curve order (tbd). Start counting at 1 to account for infinity point.
  n = 1;          # Group order (tbd). Depends on generator choice.
  MM = zeros(p);  # Curve map
  MMdelt = MM;    # Used in hiliting subgroup elements
  P0 = [];        # First non-zero point found (tbd)
  G = [];         # Generator (tbd)
  Pinf = [0 0];   # Use this to represent infinity. (May need to choose smarter if is a legit curve point)
  for x = X
    yy = f(x);
    ystr = "[";
    for y = yy
      if (length(P0)==0 && y!=0) P0 = [x,y]; end
      if (length(G)==0)
        if (GenIndex == 1)
          G = [x, y];
        else
          GenIndex--;
        end
      end
      if (length(ystr)!=1) ystr = sprintf("%s%s",ystr,", "); end
      ystr = sprintf("%s%g",ystr,y);
      N++;
      MM(y+1, x+1) = 1;   
    end
    ystr=sprintf("%s%s",ystr,"]");
    if (true || x<=30)
      printf("For x = %g, {y} = %s\n", x, ystr);
    end
  end
  printf("Field order p = %g\n", p);
  printf("Curve order N = %g\n", N);
  printf("Representing infinity as = %g, %g\n", Pinf);
  printf("Picked Generator G = %g, %g\n", G);

  P = G; PP = [G];
  printf("G *  1  =  %6s  (%s)\n", PointToString(G), Psignum(G));
  MM(P(2)+1,P(1)+1) += 1;
  MMdelt(P(2)+1,P(1)+1) += 1;      
  for i = 2:OrderGuess
    P = PointAdd(P, G);
    PP = [PP; P];
    n++;
    if(!isPinf(P))
      MM(P(2)+1,P(1)+1) += 1;
      MMdelt(P(2)+1,P(1)+1) += i;      
    end
    printf("G * %2g  =  %6s  (%s)\n", i, PointToString(P), Psignum(P));
    if(isPinf(P))
      break;
    end
  end
  MM = MM + MMdelt/n;
  printf("Params (p, a, b) = (%g, %g, %g);  G = (%g, %g)\n", p, a, b, G);
  printf("Group order:  n = %g  (cofactor N/n = h: %g/%g = %g)\n", n, N, n, N/n);
  printf("Ratio to p: n/p = %g/%g = %g\n", n, p, n/p);
  if (length(Warnings)>0)
    printf(Warnings);
  end

  Title = sprintf("Elliptic Curve with n = %g points (incl. Infinity) over F_{p} with p = %g^{ }", N, p);
  Title = sprintf("%s\nCurve Equation:  y^2 = x^3 + ax + b;  a = %g, b = %g", Title, a, b);
  XLabel = sprintf("Generator point G = (%g, %g)", G);

  imagesc(X, X, MM);
  #colorbar("eastoutside");
  axis("xy");
  title(Title);
  xlabel(XLabel, "fontweight", "bold");

  # Label Selected Curve Points:
  if (false)   # Point labels will require case-by-case adjustments so this code block defaults to skipped
    # Label Generator point:
    text(G(1)-.05, G(2)+.7, "G", "color", "white", "fontsize", 12, "verticalalignment", "bottom",
                                 "horizontalalignment", "center", "fontweight", "bold");
    # Label select multiples of Generator:
    iselect = [2, 3, 4, 5, 6, 7, 8, 9, length(PP(:,1))-3, length(PP(:,1))-2, length(PP(:,1))-1];
    for i = iselect
      dy = 0.8;
      valign="bottom";
      if (sum(i == [8,9,2,44,46])==1)   # Special cases (adjust positions)
        dy = -0.9;
        valign="top";
      end
      Ptmp = PP(i,:);
      text(Ptmp(1)-.05, Ptmp(2)+dy, sprintf("%g*G",i), "color", "white", "fontsize", 12, "verticalalignment", valign,
           "horizontalalignment", "center", "fontweight", "bold");
    end
  end
  
end


function st = PointToString(P)
  st = "(inf)";
  if (!isPinf(P))
    st = sprintf("%2g, %2g", P);
  end
end

function b = isPinf(P)
  global Pinf;
  b = false;
  if (sum(P==Pinf)==2)
    b = true;
  end
end

function s = Psignum(P)
  global param;
  p = param(1);
  symb = "-0+";
  p2 = floor(p/2);
  y = mod(P(2)+p2,p)-p2;
  s = symb(sign(y)+2);
end

# Returns empty [], [0], or [y1, y2] solutions to {y} = f(x)
function y = f(x)
  global param;
  global Xp2;   # Will use as a lookup table

  p = param(1);
  a = param(2);
  b = param(3);

  y2 = mod(x^3 + a*x + b, p); # y^2
  y = find(Xp2==y2)-1;        # Finds the sqroots of y2

end

# Point addition, Wikipedia algorithm
function Q = PointAdd(P1,P2)
  global param;
  global Pinf;
  p = param(1);
  a = param(2);

  if (isPinf(P1))
    Q = P2;
    return
  end
  if (isPinf(P2))
    Q = P1;
    return
  end

  lamnum = mod(P2(2) - P1(2),p);
  lamden = mod(P2(1) - P1(1),p);
  lamdeninv = xinv(lamden); # can be []
  if (lamnum==0 && lamden == 0) # spcl case of point doubling
    lamnum = mod(3*P1(1)*P1(1) + a, p);
    lamden = mod(2*P1(2), p);
    lamdeninv = xinv(lamden); # can be still be []
  end
  if (length(lamdeninv)==0) # Not 100% sure this is correct criteria
    Q = Pinf;
    return
  end
  lam = mod(lamnum*lamdeninv, p);
  lam2 = mod(lam*lam,p);
  Qx = mod(lam2-P1(1)-P2(1), p);
  Qy = mod(lam*(P1(1)-Qx)-P1(2), p);
  Q = [Qx, Qy];

end

# Returns inverse in Fp, else [] if inverse doesn't exist
# I.e., find xi s.t. xi * x mod p = 1
function xi = xinv(x)
  global param;
  global X;
  p = param(1);
  xi = find(mod(x*X,p)==1)-1;
end
