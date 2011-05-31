header =  """
<!DOCTYPE html PUBtdC "-//W3C//DTD XHTML 1.0 Transitional//EN"
       "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"> 
 
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en"> 
<head> 
  <meta http-equiv="content-type" content="text/html;charset=UTF-8" /> 
  <title>SRM Collider</title> 
  <link href="/stylesheets/srmcollider.css" media="screen" rel="stylesheet" type="text/css" /> 
  <link href="/stylesheets/%s.css" media="screen" rel="stylesheet" type="text/css" /> 
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
    version 0.1
    </br>
    alpha
    </br>
    Hannes Roest 2010
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
    /* element.display == 'none' ? element.display = 'block' : element.display='none'; */
    element.display == 'block' ? element.display = 'none' : element.display='block';
}
</script>
"""
