function R = JAKSTATcrossEndo_rates(x,par,t,stm)
R=[par(1)*x(32)*x(1);
par(2)*x(2)*x(2);
par(3)*x(3);
par(4)*x(4);
par(4)*x(5);
par(4)*x(6);
par(4)*x(7);
par(4)*x(8);
par(4)*x(9);
par(4)*x(10);
par(4)*x(11);
par(4)*x(12);
par(4)*x(13);
par(9)*x(34)*x(16);
par(11)*x(34)*x(1);
par(8)*x(2)*x(17);
par(3)*x(18);
par(4)*x(19);
par(4)*x(20);
par(4)*x(21);
par(4)*x(22);
par(4)*x(23);
par(4)*x(24);
par(4)*x(25);
par(4)*x(26);
par(4)*x(27);
par(4)*x(28);
par(13)*stm(1)*x(31);
par(14)*x(32);
par(15)*x(31);
par(16);
par(17)*stm(2)*x(33);
par(18)*x(34);
par(19)*x(33);
par(20);
];
end

