header =  """
<!DOCTYPE html PUBtdC "-//W3C//DTD XHTML 1.0 Transitional//EN"
       "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"> 
 
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en"> 
<head> 
  <meta http-equiv="content-type" content="text/html;charset=UTF-8" /> 
  <title>SRM Collider</title> 
  <link href="/stylesheets/srmcollider.css" media="screen" rel="stylesheet" type="text/css" /> 
  <link href="/stylesheets/%s.css" media="screen" rel="stylesheet" type="text/css" /> 

  <script src="http://www.srmcollider.org//lightbox/js/jquery-1.7.2.min.js"></script>
  <script src="http://www.srmcollider.org//lightbox/js/lightbox.js"></script>
  <link href="http://www.srmcollider.org/lightbox/css/lightbox.css" rel="stylesheet" />

</head> 

<body>
<div class="whole">
""" % 'fancy'
#options for the css are
## first
## second
## brown (quite nice)
## ocker 
## green

topdiv = """
<div class="top">
    <div class="header">
    SRM Collider
    </div>
    <div class="version">
    version 1.4
    </br>
    Hannes R&ouml;st 2012
    </div>
    <div class="navigation">
    %s
    </div>
</div>
"""


warm_welcome = topdiv % """
    <span class="active-nav"> 
        <a href="srmcollider.py">Collider </a>
    </span>
    <span class="inactive-nav"> <a href="download.html">Download </a> </span>
    <span class="inactive-nav"> <a href="about.html">About </a> </span>
"""



# Javascript function to toggle a div
toggleDisplay = """
<script language="javascript">
function toggleDisplay(e){
    element = document.getElementById(e).style;
    element.display == 'none' ? element.display = 'block' : element.display='none'; 
    /*
    element.display == 'block' ? element.display = 'none' : element.display='block';
    */
}
</script>
"""



resultInterpretation = """<p>How to interpret the results: <br>
    <ul>
    <li>
    The CSV file contains the full information with all combinations of
    transitions that are UIS (up to the selected order). 
    </li>
    <li>
    For each peptide, the website outputs a summary with 
    <ul>
    <li> The peptide sequence, q1 and SSRCalc.</li>
    <li> The number of transitions without any interference (UIS<sub>1</sub>).</li>
    <li> The peptides that interfere with the target peptide (Interfering
    peptides). Here, for each interfering peptide the transitions of the
    target peptide that are affected are listed. This means that all the
    listed combinations of transitions are <b>not</b> UIS (since they are
    shared with another peptide).</li>
    <li> The transitions of the target peptide with the interfering
    peptides listed for each transition.</li>
    </ul>
    </li>
    </ul>
    </p>"""


