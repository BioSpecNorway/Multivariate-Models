function [saisir] = Concad_Blocks(saisir,opt)
%center				- concad all blocks together

if(opt==2)
    saisir=appendcol(saisir(1),saisir(2));
elseif(opt==3)
    saisir12=appendcol(saisir(1),saisir(2));
    sasir123=appendcol(saisir12,saisir(3));
    saisir=[];
    saisir=sasir123;
elseif(opt==4)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir=[];
    saisir=appendcol(saisir12,saisir34);
elseif(opt==5)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir1234=appendcol(saisir12,saisir34);
    saisir12345=appendcol(saisir1234,saisir(5));
    saisir=[];
    saisir=saisir12345;
elseif(opt==6)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir56=appendcol(saisir(5),saisir(6));
    saisir1234=appendcol(saisir12,saisir34);
    saisir123456=appendcol(saisir1234,saisir56);
    saisir=[];
    saisir=saisir123456;
elseif(opt==7)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir56=appendcol(saisir(5),saisir(6));
    saisir1234=appendcol(saisir12,saisir34);
    saisir123456=appendcol(saisir1234,saisir56);
    saisir1234567=appendcol(saisir123456,saisir(7));
    saisir=[];
    saisir=saisir1234567;
elseif(opt==8)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir56=appendcol(saisir(5),saisir(6));
    saisir78=appendcol(saisir(7),saisir(8));
    saisir1234=appendcol(saisir12,saisir34);
    saisir123456=appendcol(saisir1234,saisir56);
    saisir12345678=appendcol(saisir123456,saisir78);
    saisir=[];
    saisir=saisir12345678;    
elseif(opt==9)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir56=appendcol(saisir(5),saisir(6));
    saisir78=appendcol(saisir(7),saisir(8));
    saisir1234=appendcol(saisir12,saisir34);
    saisir123456=appendcol(saisir1234,saisir56);
    saisir12345678=appendcol(saisir123456,saisir78);
    saisir123456789=appendcol(saisir12345678,saisir(9));
    saisir=[];
    saisir=saisir123456789;
elseif(opt==10)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir56=appendcol(saisir(5),saisir(6));
    saisir78=appendcol(saisir(7),saisir(8));
    saisir910=appendcol(saisir(9),saisir(10));
    saisir1234=appendcol(saisir12,saisir34);
    saisir123456=appendcol(saisir1234,saisir56);
    saisir12345678=appendcol(saisir123456,saisir78);
    saisir12345678910=appendcol(saisir12345678,saisir910);
    saisir=[];
    saisir=saisir12345678910;
elseif(opt==11)
    saisir12=appendcol(saisir(1),saisir(2));
    saisir34=appendcol(saisir(3),saisir(4));
    saisir56=appendcol(saisir(5),saisir(6));
    saisir78=appendcol(saisir(7),saisir(8));
    saisir910=appendcol(saisir(9),saisir(10));
    saisir1234=appendcol(saisir12,saisir34);
    saisir123456=appendcol(saisir1234,saisir56);
    saisir12345678=appendcol(saisir123456,saisir78);
    saisir12345678910=appendcol(saisir12345678,saisir910);
    saisir1234567891011=appendcol(saisir12345678910,saisir(11));
    saisir=[];
    saisir=saisir1234567891011;
end


