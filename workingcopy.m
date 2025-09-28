

function [] = main1(constellation,M)
    tic
    %generate constellation
    x = [-2:0.01:2];
    M = str2double(M);


    if constellation == "MPSK"
        z = MPSK(M);
    elseif constellation == "MQAM"
        z = MQAM(M);
       
    end
    
    
    
    %transmit symbols selected uniformly from the constellation
    
 
    if constellation == "MPSK"
        snr_db =0:0.5:25;
        [y,x] = Tx(z,snr_db,1000000);
        Pe = zeros(1,length(snr_db));
        for i = 1 : length(snr_db)
    
            
            y_hat = psk_dec(y(i,:),z,M);    
            Pe(i) = sum(y_hat ~= x)/(1000000);
            if(Pe(i) == 0)
                break;
            end
        end
        Pe_act = psk_act_err(snr_db,M);

    elseif constellation == "MQAM"
        snr_db =0:0.5:25;
        [y,x] = Tx(z,snr_db,500000);
        Pe = zeros(1,length(snr_db));
        for i = 1 : length(snr_db)
    
            
            y_hat = MLD(y(i,:),z);    
            
            Pe(i) = sum(y_hat ~= x)/(5*100000);
            if(Pe(i) == 0)
                break;
            end
        end
        Pe_act = qam_act_err(snr_db,M);
    end
        
        
    
    
    decision_region(z);
    
    figure
    semilogy(snr_db,Pe_act)
    hold on
    scatter(snr_db,Pe,'*')
    hold off
    xlabel("SNR in dB")
    ylabel("Probability")
    legend("exact error probability","Monte Carlo Error ")
    title("SEP vs Es/N0")
    
    
    
     toc
end

function [tx_symbols,tx_symbols_clean] = Tx(symb_set,snr_db,sample_size)

    symb_set = reshape(symb_set.',1,[]);
    tx_symbols_clean = datasample(symb_set,sample_size);
    
    %add complex gaussian noise
    snr = 10.^(snr_db/10);
    n = randn(1,sample_size)+ randn(1,sample_size)*1i;
    sigma = sqrt( 1./(2*snr));
    tx_symbols = tx_symbols_clean + sigma' * n;
    %scatter(real(tx_symbols),imag(tx_symbols))
    
end



function constellation = MPSK(M)
    
    k = 1 : (M);
    constellation = exp(2i*pi*(k-1)/M);
    
     
end

function constellation = MQAM(M)
    if log2(M) ~= floor(log2(M))
        constellation = -1;
        disp("M is not a power of 2")
        
    else
        r = sqrt(M);

        if r ~= floor(r)
            r = sqrt(2*M);
        end
       
        c = M/r;
        x = [0 : r-1]*2 - (r-1);
        y = [0 : c-1]*2 - (c-1);

        constellation =   x + 1i*y';      
        constellation = constellation * (1/norm(constellation,2))* sqrt(M);
        constellation = reshape(constellation.',1,[]);
      
         
        
    end
 
end

function y_hat = MLD(tx_symbols,constellation)
    
    %y_hat = zeros(1,3*100000);
    %temp = zeros(M,3*100000);
    
    tx_symbols_real = real(tx_symbols);
    tx_symbols_imag = imag(tx_symbols);

    constellation_real = real(constellation);
    constellation_imag = imag(constellation);
    %disp(size(tx_symbols_real))
    %disp(size(constellation_real))
    temp_real = (tx_symbols_real - constellation_real');
    
    temp_imag = (tx_symbols_imag - constellation_imag');
    dis = temp_real.^2 + temp_imag.^2;
    [~,k]=min(dis);
    y_hat=constellation(k);
    
    %disp(size(temp1))
    
    
    %disp(y_hat(5))
end


function [y_hat] = psk_dec(y,constellation,M)
ang = angle(y);
    neg = ang <0;
    ang(neg) = ang(neg) + 2*pi;
    ang = ang +(pi/M);
    i = (floor(ang*M/(2*pi)));
    
    y_hat = constellation(mod(i,M)+1);

end

function [y_hat] = qam_dec(y,constellation,x)
    c_real = real(constellation);
    c_imag =imag(constellation);
    y_real = real(y);
    y_imag =imag(y);

    min_real = min(c_real);
    min_imag = min(c_imag);
    
    c_real = c_real + -min_real;
    c_imag = c_imag + -min_imag;
    y_real = y_real + -min_real;
    y_imag = y_imag + -min_imag;
    
    max_real = max(c_real);
    max_imag = max(c_imag);
    
    c_real = c_real/max_real;
    c_imag = c_imag/max_imag;
    y_real = y_real/max_real;
    y_imag = y_imag/max_imag;

    r = sum(c_real ==0);
    c = sum(c_imag ==0);

    c_real = c_real*(c-1);
    c_imag = c_imag*(r-1);
    y_real = y_real*(c-1);
    y_imag = y_imag*(r-1);

    y_real = y_real+0.5;
    y_imag = y_imag+0.5;
    
    for i =1 : length(y_real)
        if y_real(i) <= 0
            y_real(i) = 0;
        elseif y_real(i) >= c-1
            y_real(i) = c-1;
        else
            y_real(i) = floor(y_real(i));
        end
    
        if y_imag(i) <= 0
            y_imag(i) = 0;
        elseif y_imag(i) >= r-1
            y_imag(i) = r-1;
        else
            y_imag(i) = floor(y_imag(i));
        end
    end

    y_real = c_real/(c-1);
    y_imag = c_imag/(r-1);

    y_real = c_real*max_real;
    y_imag = c_imag*max_imag;

    y_real = c_real + min_real;
    y_imag = c_imag + min_imag;

    y_hat = NaN;
    scatter(y_real,y_imag,"*")
    hold on
    scatter(real(x),imag(x))
    hold off

end

function theoreticalSER= psk_act_err(snr_db,M)
snr_value = 10.^(snr_db/10);
    for i = 1 : length(snr_db)
        if(M==2)
            theoreticalSER(i) = 0.5*(erfc(sqrt(snr_value(i))*sin(pi/M)));
        else
        theoreticalSER(i) = (erfc(sqrt(snr_value(i))*sin(pi/M)));
        end
    
    end
end




function theoritical_ser=qam_act_err(snr_db,M)
snr_value = 10.^(snr_db/10);

    for i = 1 : length(snr_db)
        if(M==8)
           theoritical_ser(i) = (2.5)*qfunc(sqrt(snr_value(i)/3));

        
        else
        
            theoritical_ser(i) = 2*(1-1/sqrt(M)).*erfc(sqrt(3/2.*snr_value(i)/(M-1))) - ((1-(2/sqrt(M))+1/M).*erfc(sqrt(3/2.*snr_value(i)/(M-1))).*erfc(sqrt(3/2.*snr_value(i)/(M-1))));
        
        
        end
    end
end


function [gridx] = decision_region(constellation)
    
    x = [-2:0.01:2];
    [gridy,gridx] = meshgrid(x,x);
    class = gridx;
    %temp_real = gridx;
    %temp_imag = gridx;

    constellation_real = real(constellation);
    constellation_imag = imag(constellation);
    for i = 1 : length(x)
        temp_real = (gridx(i,:) - constellation_real');
        
        temp_imag = (gridy(i,:) - constellation_imag');
        
        dis = temp_real.^2 + temp_imag.^2;
        
        [~,class(i,:)] = min(dis);
        
        %display(dis)

        %y_hat=constellation(k);
    end
    %imagesc(x,x,class)
    gridx = reshape(gridx.',1,[]);
    gridy = reshape(gridy.',1,[]);
    class = reshape(class.',1,[]);

    figure
    scatter(gridx,gridy,[],class)
    hold on
    scatter(real(constellation),imag(constellation),'filled');
    hold off
    xlabel("Real Axis")
    ylabel("Imaginary Axis")
    title("Decision Region")
    %disp(gridx)
end