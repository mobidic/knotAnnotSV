
$(document).ready(function () {

     $('#tab').DataTable();

});

/*click to change tab METADATA, ALL,etc  */
function openCity(evt, cityName) {
	var i, tabcontent, tablinks;
	tabcontent = document.getElementsByClassName("tabcontent");
	for (i = 0; i < tabcontent.length; i++) {
		tabcontent[i].style.display = "none";
	}
	tablinks = document.getElementsByClassName("tablinks");
	for (i = 0; i < tablinks.length; i++) {
		tablinks[i].className = tablinks[i].className.replace(" active", "");
	}
	document.getElementById(cityName).style.display = "block";
	evt.currentTarget.className += " active";
}


/*click to display comment*/
function myFunction() {
	  var x = document.getElementById("contentA");
	    if (x.style.display === "none") {
			    x.style.display = "block";
				  } else {
					      x.style.display = "none";
						    }
}
