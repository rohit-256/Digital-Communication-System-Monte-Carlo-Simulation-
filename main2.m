function main2(p)
tic

read_Symbol_constellation = table2array(readtable(p));
snr_db = 0:1:20;
symbolpoints = read_Symbol_constellation(:,1) + 1i*read_Symbol_constellation(:,2);
symbolpoints = symbolpoints';
[tx_symbols,x]= Tx(symbolpoints,snr_db);
 size(tx_symbols);
 
 Pe = zeros(1,length(snr_db));
 for i = 1 : length(snr_db)
    y_hat = MLD(tx_symbols(i,:),symbolpoints);

    Pe(i) = sum(y_hat ~= x)/(1000000);
 end

%---plotting SEP  vs Es/No------%
figure
semilogy(snr_db,Pe)
xlabel('SNR in dB')
ylabel('PROBABILITY')
title('SEP vs Es/No')

%Plotting Decision Region
decision_region(symbolpoints);



toc
end

%uniformly sample from the constellation and add noise
function [tx_symbols,tx_symbols_clean] = Tx(symb_set,snr_db)

    symb_set = reshape(symb_set.',1,[]);
    tx_symbols_clean = datasample(symb_set,1000000);
    
    snr = 10.^(snr_db/10);
    n = randn(1,1000000)+ randn(1,1000000)*1i;
    sigma = sqrt( 1./(2*snr));
    tx_symbols = tx_symbols_clean + sigma' * n;
    
end


function y_hat = MLD(tx_symbols,constellation)
    
    
    
    tx_symbols_real = real(tx_symbols);
    tx_symbols_imag = imag(tx_symbols);

    constellation_real = real(constellation);
    constellation_imag = imag(constellation);
    
    temp_real = tx_symbols_real - constellation_real';
    
    temp_imag = tx_symbols_imag - constellation_imag';
    dis = temp_real.^2 + temp_imag.^2;
    [~,k]=min(dis);
    y_hat=constellation(k);
    
    
end


function [gridx] = decision_region(constellation)
   
    x = [-4:0.01:4];
    [gridy,gridx] = meshgrid(x,x);
    class = gridx;
    
    constellation_real = real(constellation);
    constellation_imag = imag(constellation);
    for i = 1 : length(x)
        temp_real = (gridx(i,:) - constellation_real');
       
        temp_imag = (gridy(i,:) - constellation_imag');
       
        dis = temp_real.^2 + temp_imag.^2;
       
        [~,class(i,:)] = min(dis);
       
        
    end

    gridx = reshape(gridx.',1,[]);
    gridy = reshape(gridy.',1,[]);
    class = reshape(class.',1,[]);
    
    figure
    scatter(gridx,gridy,[],class)
    hold on
    scatter(real(constellation),imag(constellation),'filled')
    hold off
    xlabel('real axis');
ylabel('imaginary axis');
title('Decision region')
end