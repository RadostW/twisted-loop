#include"../libs/writhe.h"
#include"math.h"

namespace wr
{

wr::Curve getInitialCondition(double link)
{
    auto c = wr::Curve(15);

    bool HighWrInit = false;
    bool LowWrInit = false;

    if(HighWrInit) // High writhe init
    {
        auto c = wr::Curve(20);
        c[ 0  ]=wr::Point( 0,0,0.005674587872752547 );
        c[ 1  ]=wr::Point( 0.05422472243718961,0.0033727394280070825,-0.0038854522678721626 );
        c[ 2  ]=wr::Point( 0.11098927081720432,-0.004036227131297148,0.0018301560527810888 );
        c[ 3  ]=wr::Point( 0.17165466845698277,0.021474546564462765,0.01216848695331357 );
        c[ 4  ]=wr::Point( 0.21754814250044635,0.02430467114005569,0.00942530273548041 );
        c[ 5  ]=wr::Point( 0.23422111356277375,0,0 );
        c[ 6  ]=wr::Point( 0.21754814250044635,-0.02430467114005569,-0.00942530273548041 );
        c[ 7  ]=wr::Point( 0.17165466845698277,-0.021474546564462765,-0.01216848695331357 );
        c[ 8  ]=wr::Point( 0.11098927081720432,0.004036227131297148,-0.0018301560527810888 );
        c[ 9  ]=wr::Point( 0.05422472243718961,-0.0033727394280070825,0.0038854522678721626 );
        c[10  ]=wr::Point( 0,0,-0.005674587872752547 );
        c[11  ]=wr::Point( -0.05422472243718961,0.0033727394280070825,0.0038854522678721626 );
        c[12  ]=wr::Point( -0.11098927081720432,-0.004036227131297148,-0.0018301560527810888 );
        c[13  ]=wr::Point( -0.17165466845698277,0.021474546564462765,-0.01216848695331357 );
        c[14  ]=wr::Point( -0.21754814250044635,0.02430467114005569,-0.00942530273548041 );
        c[15  ]=wr::Point( -0.23422111356277375,0,0 );
        c[16  ]=wr::Point( -0.21754814250044635,-0.02430467114005569,0.00942530273548041 );
        c[17  ]=wr::Point( -0.17165466845698277,-0.021474546564462765,0.01216848695331357 );
        c[18  ]=wr::Point( -0.11098927081720432,0.004036227131297148,0.0018301560527810888 );
        c[19  ]=wr::Point( -0.05422472243718961,-0.0033727394280070825,-0.0038854522678721626 );
        return c;
    }
    else if(LowWrInit || link < 1.875)
    {
        //Low writhe init
        c[0]=wr::Point(0.,0.,0.1820694908908002);
        c[1]=wr::Point(0.07249260641999708,0.00809101247729739,0.16376108676623502);
        c[2]=wr::Point(0.12471906063910566,0.014194820862509904,0.12189575514685805);
        c[3]=wr::Point(0.1546812627669808,0.01793830762115063,0.06517443233674541);
        c[4]=wr::Point(0.16142893041893158,0.01878214539367789,0.002883383685381987);
        c[5]=wr::Point(0.1422637490305791,0.017637930023908666,-0.062365392387916024);
        c[6]=wr::Point(0.09972819066434507,0.012764976681821964,-0.11349831641225994);
        c[7]=wr::Point(0.03554356069257335,0.0046223749916138365,-0.14455266618872875);
        c[8]=wr::Point(-0.03554356069257335,-0.0046223749916138365,-0.14455266618872875);
        c[9]=wr::Point(-0.09972819066434507,-0.012764976681821964,-0.11349831641225994);
        c[10]=wr::Point(-0.1422637490305791,-0.017637930023908666,-0.062365392387916024);
        c[11]=wr::Point(-0.16142893041893158,-0.01878214539367789,0.002883383685381987);
        c[12]=wr::Point(-0.1546812627669808,-0.01793830762115063,0.06517443233674541);
        c[13]=wr::Point(-0.12471906063910566,-0.014194820862509904,0.12189575514685805);
        c[14]=wr::Point(-0.07249260641999708,-0.00809101247729739,0.16376108676623502);
    }
    else if(link < 2.125)
    {
        c[0 ] =wr::Point(0.,0.,-0.004936819003746011);
        c[1 ] =wr::Point(0.08469426624177992,0.06961264985843674,-0.02331468187622865);
        c[2 ] =wr::Point(0.15653534793393722,0.08239056128513103,-0.03432094556065397);
        c[3 ] =wr::Point(0.2029544841577551,0.0523321798441953,-0.028890112231375102);
        c[4 ] =wr::Point(0.21096097139130443,0.008845538068638113,-0.017950106492054887);
        c[5 ] =wr::Point(0.19111273146207403,-0.027213586890279574,-0.009334052485369452);
        c[6 ] =wr::Point(0.13949448276267273,-0.046590851053602375,-0.006702194435593414);
        c[7 ] =wr::Point(0.04568552568145135,-0.024775895130852113,-0.014842129892107831);
        c[8 ] =wr::Point(-0.04568552568145135,0.024775895130852113,-0.014842129892107831);
        c[9 ] =wr::Point(-0.13949448276267273,0.046590851053602375,-0.006702194435593414);
        c[10 ] =wr::Point(-0.19111273146207403,0.027213586890279574,-0.009334052485369452);
        c[11 ] =wr::Point(-0.21096097139130443,-0.008845538068638113,-0.017950106492054887);
        c[12 ] =wr::Point(-0.2029544841577551,-0.0523321798441953,-0.028890112231375102);
        c[13 ] =wr::Point(-0.15653534793393722,-0.08239056128513103,-0.03432094556065397);
        c[14 ] =wr::Point(-0.08469426624177992,-0.06961264985843674,-0.02331468187622865);
    }
    else if(link < 2.375)
    {
        c[0 ] =wr::Point(0.,0.,0.010597521378430065);
        c[1 ] =wr::Point(0.07859865955796944,0.061555557815427805,-0.010898915599812167);
        c[2 ] =wr::Point(0.14255604101301217,0.07874753088800308,-0.02356842482236935);
        c[3 ] =wr::Point(0.19097418073622113,0.062362207719930965,-0.01859906958391496);
        c[4 ] =wr::Point(0.21426955441643383,0.024802981077257563,-0.0041956434394351955);
        c[5 ] =wr::Point(0.1998982991128268,-0.018574200689984603,0.010485696870629184);
        c[6 ] =wr::Point(0.14396602667057018,-0.044547033405818406,0.016451838809701372);
        c[7 ] =wr::Point(0.04732993913565975,-0.024547168794872995,0.004379965898078732);
        c[8 ] =wr::Point(-0.04732993913565975,0.024547168794872995,0.004379965898078732);
        c[9 ] =wr::Point(-0.14396602667057018,0.044547033405818406,0.016451838809701372);
        c[10 ] =wr::Point(-0.1998982991128268,0.018574200689984603,0.010485696870629184);
        c[11 ] =wr::Point(-0.21426955441643383,-0.024802981077257563,-0.0041956434394351955);
        c[12 ] =wr::Point(-0.19097418073622113,-0.062362207719930965,-0.01859906958391496);
        c[13 ] =wr::Point(-0.14255604101301217,-0.07874753088800308,-0.02356842482236935);
        c[14 ] =wr::Point(-0.07859865955796944,-0.061555557815427805,-0.010898915599812167);
    }
    else if(link < 2.625)
    {
        c[0 ] =wr::Point(0.,0.,0.002873046809683538);
        c[1 ] =wr::Point(0.10622210409950755,0.02646506969074701,-0.039210553260539696);
        c[2 ] =wr::Point(0.17430456517406884,0.035419671178696426,-0.04650107876826464);
        c[3 ] =wr::Point(0.2062819835074886,0.0335823297821948,-0.020103659185154894);
        c[4 ] =wr::Point(0.20975552318219104,0.02650666321543957,0.01420447927738651);
        c[5 ] =wr::Point(0.1840315013235769,0.015462336954613387,0.04206087236638306);
        c[6 ] =wr::Point(0.12434815149778317,0.0033658770382950457,0.04133460431997879);
        c[7 ] =wr::Point(0.03951141242059983,-0.002561181088204436,0.001557954454276147);
        c[8 ] =wr::Point(-0.03951141242059983,0.002561181088204436,0.001557954454276147);
        c[9 ] =wr::Point(-0.12434815149778317,-0.0033658770382950457,0.04133460431997879);
        c[10 ] =wr::Point(-0.1840315013235769,-0.015462336954613387,0.04206087236638306);
        c[11 ] =wr::Point(-0.20975552318219104,-0.02650666321543957,0.01420447927738651);
        c[12 ] =wr::Point(-0.2062819835074886,-0.0335823297821948,-0.020103659185154894);
        c[13 ] =wr::Point(-0.17430456517406884,-0.035419671178696426,-0.04650107876826464);
        c[14 ] =wr::Point(-0.10622210409950755,-0.02646506969074701,-0.039210553260539696);
    }
    else if(link < 2.875)
    {
        c[0 ] =wr::Point(0.,0.,0.03503607257736719);
        c[1 ] =wr::Point(0.10718526866159869,0.026631278533702812,-0.005706505215611413);
        c[2 ] =wr::Point(0.17761789062794747,0.03768495281560691,-0.01176244255363526);
        c[3 ] =wr::Point(0.20858367294839036,0.040274981041327075,0.01906300214067621);
        c[4 ] =wr::Point(0.20560499259128104,0.03715641750906279,0.05477432463434627);
        c[5 ] =wr::Point(0.1760413893360367,0.029113694899817896,0.07881352027284946);
        c[6 ] =wr::Point(0.11971734952379402,0.01437788740537213,0.07414444700447605);
        c[7 ] =wr::Point(0.03957232245959401,-0.0008949040761948155,0.03309021133927969);
        c[8 ] =wr::Point(-0.03957232245959401,0.0008949040761948155,0.03309021133927969);
        c[9 ] =wr::Point(-0.11971734952379402,-0.01437788740537213,0.07414444700447605);
        c[10] =wr::Point(-0.1760413893360367,-0.029113694899817896,0.07881352027284946);
        c[11] =wr::Point(-0.20560499259128104,-0.03715641750906279,0.05477432463434627);
        c[12] =wr::Point(-0.20858367294839036,-0.040274981041327075,0.01906300214067621);
        c[13] =wr::Point(-0.17761789062794747,-0.03768495281560691,-0.01176244255363526);
        c[14] =wr::Point(-0.10718526866159869,-0.026631278533702812,-0.005706505215611413);
    }
    else if(link < 3.125)
    {
        c[0 ] =wr::Point(0.,0.,0.024260788157493076);
        c[1 ] =wr::Point(0.10792379491160616,-0.0009312844268642572,-0.013108370806999536);
        c[2 ] =wr::Point(0.1763506764939207,-0.004412947874448749,-0.022516547618648904);
        c[3 ] =wr::Point(0.2094358101352536,-0.0032707587663200364,0.00088110286426023);
        c[4 ] =wr::Point(0.2149388088982483,-0.0004489083610568742,0.03304840953569869);
        c[5 ] =wr::Point(0.19221814491495964,0.0010570102451687635,0.06305272005590458);
        c[6 ] =wr::Point(0.13069146192619863,-0.003464299902705042,0.06591739155396345);
        c[7 ] =wr::Point(0.040200727075768244,-0.010125629792539111,0.02261868524982666);
        c[8 ] =wr::Point(-0.040200727075768244,0.010125629792539111,0.02261868524982666);
        c[9 ] =wr::Point(-0.13069146192619863,0.003464299902705042,0.06591739155396345);
        c[10 ] =wr::Point(-0.19221814491495964,-0.0010570102451687635,0.06305272005590458);
        c[11 ] =wr::Point(-0.2149388088982483,0.0004489083610568742,0.03304840953569869);
        c[12 ] =wr::Point(-0.2094358101352536,0.0032707587663200364,0.00088110286426023);
        c[13 ] =wr::Point(-0.1763506764939207,0.004412947874448749,-0.022516547618648904);
        c[14 ] =wr::Point(-0.10792379491160616,0.0009312844268642572,-0.013108370806999536);
    }
    else if(link < 3.375)
    {
        c[0 ] =wr::Point(0.,0.,-0.011004676294259446);
        c[1 ] =wr::Point(0.061272925455045354,-0.00965132515128006,-0.022570567469089666);
        c[2 ] =wr::Point(0.13708314708099764,-0.05246285352537647,-0.04283570450578866);
        c[3 ] =wr::Point(0.1940560494520588,-0.05864786729123681,-0.03262950108952487);
        c[4 ] =wr::Point(0.21505356566889378,-0.03727132360209434,-0.007942739999323774);
        c[5 ] =wr::Point(0.20029977671639593,-0.00946322891763287,0.015202062397172297);
        c[6 ] =wr::Point(0.14060374091060024,0.0030106102108704647,0.018637356521480658);
        c[7 ] =wr::Point(0.04251560729940221,-0.010747121970385366,-0.012982233562005318);
        c[8 ] =wr::Point(-0.04251560729940221,0.010747121970385366,-0.012982233562005318);
        c[9 ] =wr::Point(-0.14060374091060024,-0.0030106102108704647,0.018637356521480658);
        c[10 ] =wr::Point(-0.20029977671639593,0.00946322891763287,0.015202062397172297);
        c[11 ] =wr::Point(-0.21505356566889378,0.03727132360209434,-0.007942739999323774);
        c[12 ] =wr::Point(-0.1940560494520588,0.05864786729123681,-0.03262950108952487);
        c[13 ] =wr::Point(-0.13708314708099764,0.05246285352537647,-0.04283570450578866);
        c[14 ] =wr::Point(-0.061272925455045354,0.00965132515128006,-0.022570567469089666);
    }
    else
    {
        c[0 ] =wr::Point(0.,0.,0.0010725273338000163);
        c[1 ] =wr::Point(0.061930112021985145,-0.002259842973273368,-0.013821075084170721);
        c[2 ] =wr::Point(0.13985805191587064,-0.03625799457333555,-0.04018572717047937);
        c[3 ] =wr::Point(0.20015012937301552,-0.030742840413885086,-0.039278047728769395);
        c[4 ] =wr::Point(0.21698754612823123,0.0008151087681632527,-0.016668120029337593);
        c[5 ] =wr::Point(0.19144007024177045,0.026597942843167154,0.00612754852172062);
        c[6 ] =wr::Point(0.12989637120876668,0.02237110133404431,0.012340138779504947);
        c[7 ] =wr::Point(0.04243655841530382,-0.005435831342769232,-0.0038323779847446917);
        c[8 ] =wr::Point(-0.04243655841530382,0.005435831342769232,-0.0038323779847446917);
        c[9 ] =wr::Point(-0.12989637120876668,-0.02237110133404431,0.012340138779504947);
        c[10 ] =wr::Point(-0.19144007024177045,-0.026597942843167154,0.00612754852172062);
        c[11 ] =wr::Point(-0.21698754612823123,-0.0008151087681632527,-0.016668120029337593);
        c[12 ] =wr::Point(-0.20015012937301552,0.030742840413885086,-0.039278047728769395);
        c[13 ] =wr::Point(-0.13985805191587064,0.03625799457333555,-0.04018572717047937);
        c[14 ] =wr::Point(-0.061930112021985145,0.002259842973273368,-0.013821075084170721);
    }

    return c;
}
}
