function v = BsplineWeigth(Deg,index,x)

switch(Deg)

    case 0
        w = x-index;
        v = Bspline(0,w);
    case 1
        v(2) = x-index(1);
        v(1) = 1-v(2);
        
    case 2
        w = x - index(2);
        v(2) = 3.0 / 4.0 - w * w;
		v(3) = (1.0 / 2.0) * (w - v(2) + 1.0);
        v(1) = 1.0 - v(3) - v(2);
        
    case 3
        w = x - index(2);
        v(4) = (1.0 / 6.0) * w * w * w;
        v(1) = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - v(4);
        v(3) = w + v(1) - 2.0 * v(2);
        v(2) = 1.0 - v(1) - v(3) - v(4);
        
    case 4
        w = x - index(3);
        w2 = w * w;
        t = (1.0 / 6.0) * w2;
        v(1) = 1.0 / 2.0 - w;
        v(1) = v(1)*v(1);
        v(1) = v(1)*(1.0 / 24.0) * v(1);
        t0 = w * (t - 11.0 / 24.0);
        t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
        v(2) = t1 + t0;
        v(4) = t1 - t0;
        v(5) = v(1) + t0 + (1.0 / 2.0) * w;
        v(3) = 1.0 - v(1) - v(2) - v(4) - v(5);
    case 5
        w = x - index(3);
        w2 = w * w;
        v(6) = (1.0 / 120.0) * w * w2 * w2;
        w2 = w2-w;
        w4 = w2 * w2;
        w = w- 1.0 / 2.0;
        t = w2 * (w2 - 3.0);
        v(1) = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - v(6);
        t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
        t1 = (-1.0 / 12.0) * w * (t + 4.0);
        v(3) = t0 + t1;
        v(4) = t0 - t1;
        t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
        t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
        v(2) = t0 + t1;
        v(5) = t0 - t1;
    case 6
        w = x - index(4);
        v(1) = 1.0 / 2.0 - w;
        v(1) = v(1) *v(1) * v(1);
        v(1) = v(1) * v(1) / 720.0;
        v(2) = (361.0 / 192.0 - w * (59.0 / 8.0 + w * (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)* (1.0 / 2.0 + w))))) / 120.0;
        v(3) = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w* (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w	* (-1.0 + w)))))) / 48.0;
        w2 = w * w;
        v(4) = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2* (21.0 / 4.0 - w2))) / 36.0;
        v(5) = (10543.0 / 960.0 + w * (289.0 / 16.0 + w	* (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w* (1.0 + w)))))) / 48.0;
        v(7) = 1.0 / 2.0 + w;
        v(7) = v(7) *v(7) * v(7);
        v(7) = v(7) *v(7) / 720.0;
        v(6) = 1.0 - v(1) - v(2) - v(3) - v(4)- v(5) - v(7);
    case 7
        w = x - index(4);
        v(1) = 1.0 - w;
        v(1) = v(1) *v(1);
        v(1) = v(1) *v(1) * v(1);
        v(1) = v(1) *(1.0 - w) / 5040.0;
        w2 = w * w;
        v(2) = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w* (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
        v(3) = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w	* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w	* (-5.0 + w))))))) / 240.0;
        v(4) = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2	* (-4.0 + w)))) / 144.0;
        v(5) = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w	* (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
        v(6) = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
        v(8) = w2;
        v(8) = v(8) *v(8) * v(8);
        v(8) = v(8) *w / 5040.0;
        v(7) = 1.0 - v(1) - v(2) - v(3) - v(4)	- v(5) - v(6) - v(8);
    case 8
        w = x - index(5);
        v(1) = 1.0 / 2.0 - w;
        v(1) = v(1) *v(1);
        v(1) = v(1) *v(1);
        v(1) = v(1) *v(1) / 40320.0;
        w2 = w * w;
        v(2) = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))* (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w	* (-3.0 + w)))) / 5040.0;
        v(3) = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w	* (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w	* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
        v(4) = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w	* (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w	* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
        v(5) = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2	* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
        v(6) = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w	* (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w	* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
        v(8) = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))	* (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w	* (3.0 + w)))) / 5040.0;
        v(9) = 1.0 / 2.0 + w;
        v(9) = v(9) *v(9);
        v(9) = v(9) *v(9);
        v(9) = v(9) *v(9) / 40320.0;
        v(7) = 1.0 - v(1) - v(2) - v(3) - v(4)	- v(5) - v(6) - v(8) - v(9);
    case 9
        w = x - index(5);
        v(1) = 1.0 - w;
        v(1) = v(1) *v(1);
        v(1) = v(1) *v(1);
        v(1) = v(1) *v(1) * (1.0 - w) / 362880.0;
        v(2) = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w	* (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w	* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
        v(3) = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w	* (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w	* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
        v(4) = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w	* (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w	* (-6.0 + w))))))))) / 4320.0;
        w2 = w * w;
        v(5) = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2	* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
        v(6) = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w	* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w	* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
        v(7) = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w	* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)	* (2.0 + w))))))) / 4320.0;
        v(8) = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w	* (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w	* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
        v(10) = w2 * w2;
        v(10) = v(10)*v(10) * w / 362880.0;
        v(9) = 1.0 - v(1) - v(2) - v(3) - v(4)	- v(5) - v(6) - v(7) - v(8) - v(10);
    case -1
        w = x-index(1);
        v(1) = oMom3(w);
        w = x-index(2);
        v(2) = oMom3(w);
        w = x-index(3);
        v(3) = oMom3(w);
        w = x-index(4);
        v(4) = oMom3(w);

end


end