function R = f(t,y,p, stimulus)

R = [
     [ (2*p(6)*y(13) + p(15)*y(28) - (p(1)*(y(35) + y(36))*(y(1)/p(2))^p(3))/((y(1)/p(2))^p(3) + 1) - (p(10)*(y(32) + y(33))*(y(1)/p(11))^p(12))/((y(1)/p(11))^p(12) + 1)) ];
     [ ((p(1)*(y(35) + y(36))*(y(1)/p(2))^p(3))/((y(1)/p(2))^p(3) + 1) - p(13)*y(2)*y(17) - 2*p(4)*y(2)^2 + (p(10)*(y(32) + y(33))*(y(1)/p(11))^p(12))/((y(1)/p(11))^p(12) + 1)) ];
     [ (p(4)*y(2)^2 - p(5)*y(3)) ];
     [ (p(5)*y(3) - p(6)*y(4)) ];
     [ (p(6)*y(4) - p(6)*y(5)) ];
     [ (p(6)*y(5) - p(6)*y(6)) ];
     [ (p(6)*y(6) - p(6)*y(7)) ];
     [ (p(6)*y(7) - p(6)*y(8)) ];
     [ (p(6)*y(8) - p(6)*y(9)) ];
     [ (p(6)*y(9) - p(6)*y(10)) ];
     [ (p(6)*y(10) - p(6)*y(11)) ];
     [ (p(6)*y(11) - p(6)*y(12)) ];
     [ (p(6)*y(12) - p(6)*y(13)) ];
     [ ((p(1)*p(6)*(y(35) + y(36))*(y(1)/p(2))^p(3))/((y(1)/p(2))^p(3) + 1) - 2*p(5)*p(6)*y(3) + (p(6)*p(10)*(y(32) + y(33))*(y(1)/p(11))^p(12))/((y(1)/p(11))^p(12) + 1)) ];
     [ (2*p(6)*p(7)*y(13) - 2*p(5)*p(7)*y(3)) ];
     [ (p(15)*y(28) - (p(7)*(y(32) + y(33))*(y(16)/p(8))^p(9))/((y(16)/p(8))^p(9) + 1)) ];
     [ ((p(7)*(y(32) + y(33))*(y(16)/p(8))^p(9))/((y(16)/p(8))^p(9) + 1) - p(13)*y(2)*y(17)) ];
     [ (p(13)*y(2)*y(17) - p(14)*y(18)) ];
     [ (p(14)*y(18) - p(15)*y(19)) ];
     [ (p(15)*y(19) - p(15)*y(20)) ];
     [ (p(15)*y(20) - p(15)*y(21)) ];
     [ (p(15)*y(21) - p(15)*y(22)) ];
     [ (p(15)*y(22) - p(15)*y(23)) ];
     [ (p(15)*y(23) - p(15)*y(24)) ];
     [ (p(15)*y(24) - p(15)*y(25)) ];
     [ (p(15)*y(25) - p(15)*y(26)) ];
     [ (p(15)*y(26) - p(15)*y(27)) ];
     [ (p(15)*y(27) - p(15)*y(28)) ];
     [ ((p(6)*p(7)*(y(32) + y(33))*(y(16)/p(8))^p(9))/((y(16)/p(8))^p(9) + 1) - 2*p(6)*p(14)*y(18)) ];
     [ (2*p(7)*p(15)*y(28) - 2*p(7)*p(14)*y(18)) ];
     [ (p(22) + p(17)*y(32) - p(19)*y(31) - p(16)*y(31)*stimulus{1}(t, p(32))) ];
     [ (p(16)*y(31)*stimulus{1}(t, p(32)) - p(18)*y(32) - p(20)*y(32) - p(17)*y(32)) ];
     [ (p(18)*y(32) - p(21)*y(33)) ];
     [ (p(29) + p(24)*y(35) - p(26)*y(34) - p(23)*y(34)*stimulus{2}(t, p(33))) ];
     [ (p(23)*y(34)*stimulus{2}(t, p(33)) - p(25)*y(35) - p(27)*y(35) - p(24)*y(35)) ];
     [ (p(25)*y(35) - p(28)*y(36)) ];
];
