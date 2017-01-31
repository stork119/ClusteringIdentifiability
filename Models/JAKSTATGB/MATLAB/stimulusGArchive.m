function result = stimulusGArchive(s,a)

t=s+10;
% if t<=10
%    result =0;
% elseif t > 13
%     result = 0;
if t > 10 && t <= 15
    
time = [ 
     10
     10.1
     11
     12
     13
     14
     14.9
     15
     ];


 GLNdata =a*[   
    0.0
    1.0
    1.0
    1.0
    1.0
    1.0
    1.0
    0.0   
    ];

    result =0;% interp1(time,GLNdata,t,'pchip');
elseif t > 25 && t <= 30
        
time = sort([25:1:30 25.1 29.9]);


 GLNdata = a*[   
    0.0
    1.0
    1.0
    1.0
    1.0
    1.0
    1.0
    0.0   
    ];

    result = interp1(time,GLNdata,t,'linear');
else
    result = 0;
end

     

end

