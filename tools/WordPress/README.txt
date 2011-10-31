Wordpress plugin for automatic annotation of TOPPAS workflows uploaded to the WordPress homepage.

Installation:
1) Copy the .php and .js file to the WP plugins directory

2) Insert the code below into the HTML section of the WP-page you want to insert .toppas file links into
   Make sure that there are no empty lines, as WP will add HTML code which destroys the JavaScript (check in final page source code)


<!-- DO NOT REMOVE THIS. Otherwise the Show/Hide() button will not work -->
<script type="text/javascript">// <![CDATA[
  $(document).ready(function() 
  {
    /* toggle (view on/off) all elements with the same class as the sending button */
    $("button").click(function () {
                        var $target = $(this);
                        var $id = "." + $target.attr("class");
                        //alert($(event.target).attr("class"));
                        $($id + "_e").toggle();
                      }
                     );
    /* hide in the beginning */
    $(".toggleClass").toggle(); 
  });// document ready
// ]]></script>


3) Insert .toppas file links preceded by four "!" and suffixed by four "#", e.g.
	!!!!http://open-ms.sourceforge.net/wp-content/uploads/2011/08/iTRAQ_Quantation_and_ID.toppas####
