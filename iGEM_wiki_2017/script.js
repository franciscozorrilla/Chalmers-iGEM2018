
	//This is for the navbar on the left of the text. Might need to change the nav part, if that conflicts with the top navbar later.
	//Note that the click function doesn't work in the small preview, but works in a separate window, for whatever reason...
	$('#side-nav nav a').on('click', function(event) {
		$("#nav li a.active").removeClass("active"); //Remove any "active" class  
		$(this).addClass("active"); //Add "active" class to selected tab 
	});

	//($(this).offset().top - VARIABLE) - VARIABLE is used to switch to the next tab
	//sooner than the top of the window reaches the section in question.
	$(window).on('scroll', function() {
	$('.target').each(function() {
		var half_height = $(window).height()/2;
		if($(window).scrollTop() >= ($(this).offset().top - half_height)) {
				var id = $(this).attr('id');
				$('#side-nav nav a').removeClass('active');
				$('#side-nav nav a[href="#'+ id +'"]').addClass('active');
		}
	});
	});


	/**
	* HEADERNAV STICKY
	* This part does the "fixed navigation after scroll" functionality
	* We use the jQuery function scroll() to recalculate our variables as the 
	* page is scrolled/
	*/
	$(window).scroll(function(){
		var window_top = $(window).scrollTop();
		var div_top = $('#header-nav-anchor').offset().top;
				if (window_top > div_top) {
						$('#header_navbar nav').addClass('header-nav_stick');
				} else {
						$('#header_navbar nav').removeClass('header-nav_stick');
				}
	});

	/**
	* SIDENAV STICKY
	* This part does the "fixed navigation after scroll" functionality
	* We use the jQuery function scroll() to recalculate our variables as the 
	* page is scrolled/
	*/
	$(window).scroll(function(){
		var window_top = $(window).scrollTop() + 0; // the number should equal the margin-top value for nav.stick
		var div_top = $('#nav-anchor').offset().top;
				if (window_top > div_top) {
						$('#side-nav nav').addClass('side-nav_stick');
						sidenav_bottom = $('#nav-anchor nav').outerHeight(true); + 0; /*62*/
				} 
				else {
						$('#side-nav nav').removeClass('side-nav_stick');
				}
	});

	/**
	* Hiding the sidenav when get too close to the footer, so that it doesn't cover it
	*/
	$(window).scroll(function(){
	 var footer_top = $('#footer-top')[0].getBoundingClientRect().top;
	 var test = 1000;
	 if (footer_top < sidenav_bottom + 200) {
			 $('#side-nav nav').addClass('hidden-sidenav');
	 }
	 else {
			 $('#side-nav nav').removeClass('hidden-sidenav');
	 }
	});

	//Smooth scrolling to links. Note the value after scrollTop: $(hash).offset().top.
	//It's necessary to compensate for the top navbar
	$(document).ready(function(){
	// Add smooth scrolling to all links
	$("a").on('click', function(event) {

		// Make sure this.hash has a value before overriding default behavior
		if (this.hash !== "") {
			// Prevent default anchor click behavior
			event.preventDefault();

			// Store hash
			var hash = this.hash;

			// Using jQuery's animate() method to add smooth page scroll
			// The optional number (800) specifies the number of milliseconds it takes to scroll to the specified area
			$('html, body').animate({
				scrollTop: $(hash).offset().top - 80 //CHANGE THIS IF TOP NAVBAR SIZE IS CHANGED
			}, 800, function(){

				// Add hash (#) to URL when done scrolling (default click behavior)
				window.location.hash = hash;
			});
		} // End if
	});
	});



	/**
	* This is the function that expands profile cards
	*/
	$(document).ready(function(){
	$('.profile-container').hover(function() {
		$(this).animate({
				width: 350,
				height: 600,
				top: -20,
				left: -25
		}, "fast");
	$(this).animate().css('box-shadow', '1px 1px 7px #000');
		$(this).css({
				zIndex: 100 
		});
	}, function() {
	$(this).animate().css('box-shadow', 'none')
		$(this).animate({
				width: 300,
				height: 350,
				top: 0,
				left: 0
		}, "fast");
		$(this).animate().css('box-shadow', '1px 1px 4px #666')
		$(this).css({
				zIndex: 1
		});
	});
	});
