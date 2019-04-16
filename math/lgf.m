function y=lgf(x)
%  Vectorized and pre-calculated version by:
%    O. Firstenberg 2012, Harvard University HQOC / MIT
%
%  Based on original code by:
%    J. Pritchard Durham University 2009
%    [downloaded from: http://massey.dur.ac.uk/jdp/code.html, Jan. 2012]

%   the pre-calculation was created using:
%   lgf_precalculated=log(factorial(0:169)); sprintf('%15.9f \n',lgf_precalculated);

lgf_precalculated=[0.000000000    0.000000000    0.693147181    1.791759469    3.178053830    4.787491743    6.579251212    8.525161361   10.604602903   12.801827480   15.104412573   17.502307846   19.987214496   22.552163853   25.191221183   27.899271384   30.671860106   33.505073450   36.395445208   39.339884187   42.335616461   45.380138898   48.471181352   51.606675568   54.784729398   58.003605223   61.261701761   64.557538627   67.889743137   71.257038967   74.658236349   78.092223553   81.557959456   85.054467018   88.580827542   92.136175604   95.719694542   99.330612455  102.968198615  106.631760261  110.320639715  114.034211781  117.771881400  121.533081515  125.317271149  129.123933639  132.952575036  136.802722637  140.673923648  144.565743946  148.477766952  152.409592584  156.360836303  160.331128217  164.320112263  168.327445448  172.352797139  176.395848407  180.456291418  184.533828861  188.628173424  192.739047288  196.866181673  201.009316399  205.168199483  209.342586753  213.532241495  217.736934114  221.956441819  226.190548324  230.439043566  234.701723443  238.978389562  243.268849003  247.572914096  251.890402210  256.221135550  260.564940972  264.921649799  269.291097651  273.673124286  278.067573440  282.474292688  286.893133295  291.323950094  295.766601351  300.220948647  304.686856766  309.164193580  313.652829950  318.152639620  322.663499127  327.185287704  331.717887197  336.261181979  340.815058871  345.379407062  349.954118041  354.539085519  359.134205370  363.739375556  368.354496072  372.979468886  377.614197874  382.258588773  386.912549123  391.575988217  396.248817052  400.930948279  405.622296161  410.322776527  415.032306728  419.750805600  424.478193418  429.214391867  433.959323995  438.712914186  443.475088121  448.245772745  453.024896238  457.812387981  462.608178527  467.412199572  472.224383927  477.044665493  481.872979230  486.709261137  491.553448223  496.405478487  501.265290892  506.132825342  511.008022665  515.890824588  520.781173716  525.679013516  530.584288294  535.496943180  540.416924106  545.344177791  550.278651724  555.220294147  560.169054037  565.124881095  570.087725725  575.057539025  580.034272767  585.017879389  590.008311976  595.005524249  600.009470555  605.020105849  610.037385686  615.061266207  620.091704128  625.128656731  630.172081848  635.221937855  640.278183660  645.340778693  650.409682896  655.484856711  660.566261076  665.653857411  670.747607612  675.847474040  680.953419514  686.065407302  691.183401114  696.307365094  701.437263809];
y=x;
y(x<170)=lgf_precalculated(x(x<170)+1);
Ilargex=x>=170;x=x(Ilargex);

%  Stirlings approximation ln(n!) = nln(n)-n+0.5ln(2pin)
y(Ilargex)=x.*log(x)-x+0.5*log(2*pi*x)+1./(12*x)-1./(360*x.^3)+1./(1260*x.^5)-1./(1680*x.^7)+1./(1188*x.^9);
