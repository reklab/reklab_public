function s = nlmtst(S);

            % test routine for spect
            disp(['nlmtst for spect']);
            figure (1);
            z=nlid_sim('L2');
            x=z(:,1);
           % Test NFFT
           for i=1:4,
                nFFT=64*i;
                s=spect(x, 'nFFT',nFFT);
                
                subplot (2,2,i);
                plot(s); set (gca,'xlim',[0 50],'ylim',[0 .06]);
                title (['NFFT=' num2str(nFFT)])
           end

           %% Test conifdence oevels
           s=spect(x,'nFFT',nFFT,'confidenceLevel',.95);
           disp(s)
           figure(2)
           subplot (3,1,1)
           plotConfidence(s);
           subplot (3,1,2);
           plotConfidence(s,'semilog');
           subplot (3,1,3);
           plotConfidence (s,'loglog');


            
            % Test cross spectrA
            figure(3);
            subplot (2,2,1);
            s=spect(x);
            plot(s); set (gca,'xlim',[0 50]);
            title (['nFFT=' num2str(nFFT)])
            title('Autospectrum of input');
            
            
            x=z(:,2);
            s=spect(x);
            subplot (2,2,2);
            plot(s);
            title('Autospectrum of output');
            
            
            
            s=spect(z);
            set(s,'comment','cross spectrum');
            subplot (2,2,3);
            plot(abs(s));
            title('Magnitude of cross spectrum;')
            subplot (2,2,4);
            plot(phase(s))
            title('Cross spectrum phase');
        end
        