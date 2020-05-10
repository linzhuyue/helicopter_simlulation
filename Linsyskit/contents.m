%Linear Systems Toolkit - Version 0.99.  Released in August 2004
%
%
%Decompositions of Autonomous Systems
% ssd       - continuous-time stability structural decomposition
% dssd      - discrete-time stability structural decomposition
% jcf       - Jordan canonical form
% rjd       - real Jordan decomposition
% 
%Decompositions of Unforced and Unsensed Systems
% osd       - observability structural decomposition
% obvidx    - observability index
% bdosd     - block diagonal observable structural decomposition
% csd       - controllability structural decomposition
% ctridx    - controllability index
% bdcsd     - block diagonal controllable structural decomposition
%
%Decompositions and Structural Properties of Proper Systems
% scbraw    - raw decomposition without integration chains
% scb       - decomposition of a continuous-time system
% dscb      - decomposition of a discrete-time system
% kcf       - Kronecker canonical form for system matrices
% morseidx  - Morse indices
% blkz      - blocking zeros
% invz      - invariant zero structure
% infz      - infinite zero structure
% l_invt    - left invertibility structure
% r_invt    - right invertibility structure
% normrank  - normal rank
% v_star    - weakly unobservable subspace
% v_minus   - stable weakly unobservable subspace
% v_plus    - unstable weakly unobservable subspace
% s_star    - strongly controllable subspace
% s_minus   - stable strongly controllable subspace
% s_plus    - unstable strongly controllable subspace
% r_star    - controllable weakly unobservable subspace
% n_star    - distributionally weakly unobservable subspace
% s_lambda  - geometric subspace S_{lambda}
% v_lambda  - geometric subspace V_{lambda}
%
%Operations of Vector Subspaces
% ssorder   - ordering of vector subspaces
% ssintsec  - intersection of vector subspaces
% ssadd     - addition of vector subspaces
%
%Decompositions and Properties of Descriptor Systems
% ea_ds     - decomposition of a matrix pair (E,A)
% sd_ds     - decomposition for descriptor systems
% invz_ds   - descriptor system invariant zero structure
% infz_ds   - descriptor system infinite zero structure
% l_invt_ds - descriptor system left invertibility structure
% r_invt_ds - descriptor system right invertibility structure
%
%System Factorizations
% mpfact    - continuous minimum-phase/all-pass factorization
% iofact    - continuous-time inner-outer factorization
% gcfact    - continuous generalized cascade factorization
% dmpfact   - discrete minimum-phase/all-pass factorization
% diofact   - discrete-time inner-outer factorization
%
%Structural Assignment via Sensor/Actuator Selection
% sa_sen    - structural assignment via sensor selection
% sa_act    - structural assignment via actuator selection
%
%Asymptotic Time-scale and Eigenstructure Assignment (ATEA)
% atea      - continuous-time ATEA
% gm2star   - infimum for continuous-time H2 control
% h2care    - solution to continuous-time H2 ARE
% h2state   - solution to continuous-time H2 control
% gm8star   - infimum for continuous-time H-infinity control
% h8care    - solution to continuous-time H-infinity ARE
% h8state   - solution to continuous-time H-infinity control
% addps     - solution to continuous disturbance decoupling
% datea     - discrete-time ATEA
% dare      - solution to general discrete-time ARE
% dgm2star  - infimum for discrete-time H2 control
% h2dare    - solution to discrete-time H2 ARE
% dh2state  - solution to discrete-time H2 control
% dgm8star  - infimum for discrete-time H-infinity control
% h8dare    - solution to discrete-time H-infinity ARE
% dh8state  - solution to discrete-time H-infinity control
% daddps    - solution to discrete-time disturbance decoupling
%
%Disturbance Decoupling with Static Output Feedback
% ddpcm     - solution to disturbance decoupling problem with static output feedback (DDPCM)
% rosys4ddp - irreducible reduced-order system that can be used to solve DDPCM
%
%
% There are 66 m-functions in the above Linear Systems Toolkit. Some of these
% m-functions are interactive, which require users to enter additional parameters
% when executed. Some can return results either in a symbolic or numerical form.
