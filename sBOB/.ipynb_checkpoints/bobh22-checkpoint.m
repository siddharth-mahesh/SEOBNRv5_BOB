%Plotting parameters
format long
font = 'Arial';
fontsize = 30;
linewidth = 3;
truncflag = false;
plotflag = true;
ploth = true;
plotpsi4 = true;

%Physical parameters
af = 0.6864;    %values for equal mass nonspinning merger, from https://arxiv.org/pdf/1305.5991.pdf
Mf = 0.9516;    %can use a fit to NR or a "first principles approach" for these
Mhalf = 1-(1-Mf)/2;  %mass after half of the total energy that will be emitted in GWs is gone, roughly appropriate for normalizing the peak strain

%Time alignment parameters
pct = 0.99; %fractional amplitude to be used for peak finding, i.e. all amplitude data > pct*Apeak is fit to a parabola to find the interpolated peak
tmin = -200; %time range for psi4, NOT h. psi4 will peak at t=0, h will ultimately be shifted to peak at t=tp
tmax = 200;
tp = 0; %time when you want strain peak amplitude to occur
dt = 0.1; %time resolution
t = (tmin:dt:tmax)';    %array of times, ranges from tmin to tmax in steps of dt

%reference frequencies, capital letters for orbital quantities, lowercase for GW quantities
Omisco = (-0.091933*af + 0.097593)/(af^2 - 2.4228*af + 1.4366); %my fit, good to O(0.1%) for af < 0.99
omqnm = 1.5251-1.1568*(1-af)^0.1292;  %Berti fit, from gr-qc/0512160
omqnm = omqnm/Mf;
Omqnm = omqnm/2;    %orbital frequency for null perturbation sourcing the l=m=2 QNM
Q = 0.7+1.4187*(1-af)^(-0.499);    %Berti fit, from gr-qc/0512160
tau = Q/omqnm;
Tau = 2*tau;    %orbital Lyapunov-based e-folding time for amplitude of l=m=2

%Amplitude
Ap = 0.908*(1-Mf).^0.794;  %my fit to SXS sims, good to O(1%) for all cases with NR errors that small
A = Ap*sech(t/Tau);    %psi4 amplitude, will peak at t=0

%Orbital frequency, see BOB paper for derivation
Omref = Omisco.*(1-af); %lim(t -> -Inf) Omega, generally 0 <= Omref <= Omisco, this choice empirically works best at early times
%NB: the choice of Omref has essentially no effect near, at, or after the
%psi4 peak, it only affects things ~10M before the peak or earlier, but
%since tp_strain < tp_psi4, this choice will affect the strain peak frequency by a few %
Omp = ((Omqnm^4 + Omref^4)/2)^(1/4);
Omm = ((Omqnm^4 - Omref^4)/2)^(1/4);
Om = (Omp^4 + Omm^4*tanh((t-tp)/Tau)).^(1/4);
%NB: if Omref = 0, then Om = Omqnm*((1 + tanh(t/Tau))/2)^(1/4);

%Orbital phase, see BOB paper for derivation
kappap = (Omp^4-Omm^4)^(1/4);
kappam = (Omp^4+Omm^4)^(1/4);
Phi = Tau*...
    ( (atan(Om/kappam) - atan(Omp/kappam))*kappam ...
    + (atanh(Om/kappam) - atanh(Omp/kappam))*kappam ...
    - (atan(Om/kappap) - atan(Omp/kappap))*kappap ...
    - (atanh(Om/kappap) - atanh(Omp/kappap))*kappap);

%GW frequency and phase from orbital values
om = 2*Om;
phi = 2*Phi;

%complex psi4, psi4 = psi4_+ + i*psi4_x
psi4 = A.*exp(1i*phi);

%complex strain, h = h_+ + i*h_x. For quasicircular motion, psi4 ~ -omega^2*h
h = -(Mhalf/Mf)^2*om.^(-2).*psi4;    %NB: if ~half the energy is radiated before the peak, this mass normalization makes the most sense
Ah = abs(h);    %assuming quasicircularity, only strain amplitude and psi4 amplitude differ, not phases
phih = unwrap(angle(h));    %at this point, phih = phipsi4, but we then shift phih to be zero at the strain peak

%calculating time, amplitude, and frequency of peak strain, shifting peak to t=tp
%the next 6 lines find the peak strain array value Aph, then fits a parabola to all
%values > pct*Aph, then uses that fit to find the interpolated "true" peak strain
% [Aph, index] = max(Ah);
% indexleft = dsearchn(Ah(1:index),pct*Aph);
% indexright = index + dsearchn(Ah(index+1:end),pct*Aph);
% [pcoeff,~,mu] = polyfit(t(indexleft:indexright),Ah(indexleft:indexright),2);
% tph = mu(1)-mu(2)*pcoeff(2)/2/pcoeff(1);    %"true" peak time from fit
% Aph = polyval(pcoeff,tph,[],mu);    %"true" peak amplitude
tph = tp + Tau*log(Omref/Omqnm);    %"true" peak time
Aph = spline(t,Ah,tph); %"true" peak amplitude
phih = phih - spline(t,phih,tph);   %making the phase=0 at the interpolated peak time
h = Ah.*exp(1i*phih);   %this strain has phase=0 at the strain peak
t = t-tph+tp;   %this shifts the time so that the peak occurs at tp, specified above

%truncating all data before the peak. Note that t(1) will NOT be exactly zero,
%because the "true" peak occurred between data points and was interpolated
if truncflag
    icut = dsearchn(t,tp)-1;
    t(1:icut) = [];
    h(1:icut) = [];
end

if plotflag
    if ploth
        figure
        plot(t,real(h),'LineWidth',linewidth)
        hold all
        plot(t,imag(h),'LineWidth',linewidth)
%         plot([tph tph],1.2*[-Aph Aph],'--k','LineWidth',linewidth)
        xlim([tmin tmax])
        xlabel('$t/M$','FontSize',fontsize,'Interpreter','latex')
        ylabel('$rh_{22}/M$','FontSize',fontsize,'Interpreter','latex')
        set(gca,'FontSize',fontsize,'LineWidth',linewidth)
        legend({'$r\Re(h_{22})/M$','$r\Im(h_{22})/M$','peak'},'FontSize',fontsize,'Interpreter','latex')
        figure
        plot(t,Ah,'LineWidth',linewidth)
        xlim([tmin tmax])
        xlabel('$t/M$','FontSize',fontsize,'Interpreter','latex')
        ylabel('$r|h_{22}|/M$','FontSize',fontsize,'Interpreter','latex')
        set(gca,'FontSize',fontsize,'LineWidth',linewidth)
    end
    if plotpsi4
        figure
        plot(t,real(psi4),'LineWidth',linewidth)
        hold all
        plot(t,imag(psi4),'LineWidth',linewidth)
        xlim([tmin tmax])
        xlabel('$t/M$','FontSize',fontsize,'Interpreter','latex')
        ylabel('$rh_{22}/M$','FontSize',fontsize,'Interpreter','latex')
        set(gca,'FontSize',fontsize,'LineWidth',linewidth)
        legend({'$rM\Re(\psi_{4,22})$','$rM\Im(\psi_{4,22})$'},'FontSize',fontsize,'Interpreter','latex')
        figure
        plot(t,A,'LineWidth',linewidth)
        xlim([tmin tmax])
        xlabel('$t/M$','FontSize',fontsize,'Interpreter','latex')
        ylabel('$rM|\psi_{4,22}|$','FontSize',fontsize,'Interpreter','latex')
        set(gca,'FontSize',fontsize,'LineWidth',linewidth)
    end
    figure
    plot(t,phi,'LineWidth',linewidth)
    xlim([tmin tmax])
    xlabel('$t/M$','FontSize',fontsize,'Interpreter','latex')
    ylabel('$\phi$','FontSize',fontsize,'Interpreter','latex')
    set(gca,'FontSize',fontsize,'LineWidth',linewidth)
    figure
    plot(t,om,'LineWidth',linewidth)
    xlim([tmin tmax])
    xlabel('$t/M$','FontSize',fontsize,'Interpreter','latex')
    ylabel('$\omega$','FontSize',fontsize,'Interpreter','latex')
    set(gca,'FontSize',fontsize,'LineWidth',linewidth)
end