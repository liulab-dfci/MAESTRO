var heatmap = null;

var giggleUrl            = "http://localhost:8080/";
var giggleUCSCBrowserUrl = "http://localhost:8081/";
var ucscBrowserUrl       = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19";

var sourceFileMap = {};	
var coordMap = {};
var def = null;
var dataForChart = null;
var valueField = "overlaps";

var uploadedFileName = null;

var ucscFileMap = {};
var ucscTrackNames = [];

$.urlParam = function(name){
    var results = new RegExp('[\?&]' + name + '=([^&#]*)').exec(window.location.href);
    if (results==null){
        return null;
    }else{
        return results[1] || 0;
    }
}

$(document).ready(function() {
        if (decodeURIComponent($.urlParam('primary_index')) != 'null') {
           giggleUrl = "http://" + decodeURIComponent($.urlParam('primary_index'));
        } 

        if (decodeURIComponent($.urlParam('ucsc_index')) != 'null') {
           giggleUCSCBrowserUrl = "http://" + decodeURIComponent($.urlParam('ucsc_index'));
        } 


	$.material.init();

	$('#giggle-url').val(giggleUrl);
	$('#giggle-tracks-url').val(giggleUCSCBrowserUrl);

	initBedUploadForm();

	promiseLoadMatrixDefinition().then( function() {
		heatmap = new heatmapD3().cellSize(15)
	                             .legendCellSize(20)
	                             .margin({top: 10, bottom: 10, left: 10, right: 70})
	                             .colors(colorbrewer.YlGnBu[9])
	                             //.colors(colorbrewer.Oranges[9]);
	                             .cellValue( function(d) { return +d[valueField]; } )
	                             .on('d3click', function(d,i) {
	                             	loadOverlapDetail(d.name, d.row, d.col);
	                             });


		loadHeatmapForRegion();

	});

        loadUCSCDefinition();
});

function loadUCSCSmartView() {

    	if ($('#overlaps').val() == null || $('#overlaps').val().trim() == "") {
		$('#no-region-warning').removeClass("hide");
		return;
	}

	var giggleTracksUrl = giggleUCSCBrowserUrl + "?region=" + $('#overlaps').val();

	$.ajax({
	    url: giggleTracksUrl,
	    type: "GET",
	    crossDomain: true,
	    dataType: "text",
	    success: function(data) {
	    	var records = [];
			data.split("\n").forEach(function(row) {
				if (row == null || row.trim() == "") {

				} else {
					fields = row.split("\t");

					
					var rec = {};
					rec.name     = fields[0].split("/")[1];
					rec.size     = fields[1];
					rec.overlaps = fields[2];

					var pos      = ucscFileMap[rec.name];
					rec.pos      = pos;
					rec.trackName = ucscTrackNames[+pos];
					records.push(rec);

				}
			});

			//var ucscTracksUrl = ucscBrowserUrl + '&position=' + chr + ":" + start + '-' + end;
			var ucscTracksUrl = ucscBrowserUrl + '&position=' + $('#overlaps').val();
			records.forEach( function(record) {
                                if (record.trackName) {
				    if (+record.overlaps > 0) {
					    ucscTracksUrl += "&" + record.trackName + "=dense";
				    } else {
					    ucscTracksUrl += "&" + record.trackName + "=hide";
                                    }
                                }
			});
			var newTab = window.open(ucscTracksUrl, '_blank');
			//newTab.focus();
	    },
	    error: function(error) {
	    	console.log("An error occurred when getting UCSC track info " + giggleTracksUrl);
	    	console.error();
	    }
	});

//function loadUCSCTracks(chr, start, end) {
}

function loadUCSCDefinition() {

	var giggleTracksDefUrl = giggleUCSCBrowserUrl + "?data";

	$.ajax({
	    url: giggleTracksDefUrl,
	    type: "GET",
	    crossDomain: true,
	    dataType: "text",
	    success: function(data) {

	    	var def = JSON.parse(data);
			def.sourceFiles.forEach( function( sourceFile ) {
				var fileName = sourceFile.name.split("/")[1];
				ucscFileMap[fileName] = sourceFile.position[0];
			});

			def.dimensions[0].elements.forEach( function(trackName) {
				ucscTrackNames.push(trackName);
			})
			

		},
	    error: function(error) {
	    	console.error;
	    }
	});	

}

function loadUCSCTracks(chr, start, end) {
	var giggleTracksUrl = giggleUCSCBrowserUrl + "?region=" + chr + ":" + start + '-' + end;
	$.ajax({
	    url: giggleTracksUrl,
	    type: "GET",
	    crossDomain: true,
	    dataType: "text",
	    success: function(data) {
	    	var records = [];
			data.split("\n").forEach(function(row) {
				if (row == null || row.trim() == "") {

				} else {
					fields = row.split("\t");
					
					var rec = {};
					rec.name     = fields[0].split("/")[1];
					rec.size     = fields[1];
					rec.overlaps = fields[2];

					var pos      = ucscFileMap[rec.name];
					rec.pos      = pos;
					rec.trackName = ucscTrackNames[+pos];
					records.push(rec);

				}
			});

			var ucscTracksUrl = ucscBrowserUrl + '&position=' + chr + ":" + start + '-' + end;
			records.forEach( function(record) {
				if (+record.overlaps > 0) {
					ucscTracksUrl += "&" + record.trackName + "=dense";
				} else {
					    ucscTracksUrl += "&" + record.trackName + "=hide";
                                }
			});
			var newTab = window.open(ucscTracksUrl, '_blank');
			//newTab.focus();




	    },
	    error: function(error) {
	    	console.log("An error occurred when getting UCSC track info " + giggleTracksUrl);
	    	console.error();
	    }
	});
}


function loadOverlapDetail(fileName, row, col) {
	var rowLabel = coordMap[row + '-' + col].rowLabel;
	var colLabel = coordMap[row + '-' + col].colLabel;

	var detailUrl = giggleUrl;
	if (uploadedFileName) {
		detailUrl += "?query=" + uploadedFileName;
	} else {
		detailUrl += "?region=" + $('#overlaps').val();
	}
	detailUrl += "&files=" + fileName + "&full";
	$.ajax({
	    url: detailUrl,
	    type: "GET",
	    crossDomain: true,
	    dataType: "text",
	    success: function(data) {
	    	var results = [];
	    	var header = null;
	    	var records = [];
			data.split("\n").forEach(function(row) {
				if (row == null || row.trim() == "") {

				} else {
					fields = row.split("\t");
					if (row.indexOf("#") == 0) {
						if (header) {
							results.push({'header': header, 'rows': records});
							recs = [];
						}
						header = {};
						header.name = fields[0].split("/")[1];
						header.size = fields[1];
						header.overlaps = fields[2];
					} else {
						var rec = {};
						rec.chr   = fields[0];
						rec.start = fields[1];
						rec.end   = fields[2];
						records.push(rec);
					}

				}
			});
			if (header) {				
				results.push({'header': header, 'rows': records});
			}
			$('#overlaps-modal .modal-header').html("<h5>" + rowLabel + "  -  " + colLabel + "</h5>");
			$('#overlaps-modal .modal-body').html("");

			results.forEach(function(result) {
				var content = "";
				content += "<table style='width:100%'>";
				
				var rowNbr = 1;
				result.rows.forEach( function(row) {
					content += 
						 "<tr>" 
						+   "<td>" + rowNbr++ + ".</td>"
						+ "<td>" + "<a href='javascript:void(0)' onclick=\"loadUCSCTracks(" + "'" + row.chr + "'," + row.start + ',' + row.end + ")\">" +row.chr + ' ' + addCommas(row.start) + '-' +  addCommas(row.end) + "</a></td>"
					    + "</tr>";
				});
				content += "</table>";
				$('#overlaps-modal .modal-body').append(content);
			});



			$('#overlaps-modal').modal('show');
			
		  

	    },
	    error: function(error) {
	    	console.log("An error occurred when getting overlap detail " + detailUrl);
	    	console.error();
	    }
	});	

}


function promiseLoadMatrixDefinition() {
	return new Promise( function(resolve, reject) {
		var defUrl  = giggleUrl + "?data";
		def = null;
		sourceFileMap = {};	

		// Get matrix definition
		$.ajax({
		    url: defUrl,
		    type: "GET",
		    crossDomain: true,
		    dataType: "text",
		    success: function(data) {
		    	def = JSON.parse(data);
				def.sourceFiles.forEach( function( sourceFile ) {
					var cellCoord = {};
					cellCoord.row = sourceFile.position[0];
					cellCoord.col = sourceFile.position[1];
					sourceFileMap[sourceFile.name] = cellCoord;
					coordMap[cellCoord.row + "-" + cellCoord.col] = {rowLabel: def.dimensions[0].elements[cellCoord.row], colLabel: def.dimensions[1].elements[cellCoord.col]};
				});
				def.cells = [];
				resolve(def);

			},
		    error: function(error) {
		    	console.error();
		    	reject('Unable to get matrix definition ' + error);
		    }
		});

	});
}

function initBedUploadForm() {


	$('#bed-upload-form').submit( function(e){
		e.preventDefault();

		$(".input-panel").removeClass("selected");

		// Validate that file was uploaded
		$('.alert').addClass("hide");
		if ($('#bed-upload-form input[type=file]')[0].files.length == 0) {
			$('#no-file-warning').removeClass("hide");		
			return;	
		}

		uploadedFileName =  $('#bed-upload-form input[type=file]')[0].files[0].name;

		// Enable odds ratio and default as checked
		//$("#radio-value-ratio").prop("checked", true);
		//$("#radio-value-ratio-span").removeClass("hide");
		$("#radio-value-combo").prop("checked", true);
		$("#radio-value-combo-span").removeClass("hide");

		// Highlight upload panel
		$("#upload-panel").addClass("selected");

		// Show loading gif
		$('.loader').removeClass("hide");


	    
	    var formData = new FormData(this);
		//formData.append(uploadedFileName, $('#bed-upload-form input[type=file]')[0].files[0]);
	    
	    getGiggleUrls();

	    var url = giggleUrl + "filepost";

	    $.ajax({
	        url         : url,
	        data        : formData,
 	        cache       : false,
	        contentType : false,
	        processData : false,
	        type        : 'POST',
	        success     : function(data, textStatus, jqXHR){
	            $('.loader').addClass("hide");
	            loadHeatmapChart(data, def);
	        },
	        error       : function(error) {
                    $('.loader').addClass("hide");
                    if (error.success().hasOwnProperty("responseText") && error.success().responseText.length > 0) {
                        loadHeatmapChart( error.success().responseText, def);
                    }
                    console.error();      	
	        }
	    });
	});	
}

function getGiggleUrls() {
	giggleUrl = $('#giggle-url').val();
	giggleUCSCBrowserUrl = $('#giggle-tracks-url').val();

	giggleUrl += '/'
	giggleUCSCBrowserUrl += '/'
}


function loadHeatmapForRegion() {
	$(".input-panel").removeClass("selected");

	// Validate that region was filled in
	$('.alert').addClass("hide");
	if ($('#overlaps').val() == null || $('#overlaps').val().trim() == "") {
		$('#no-region-warning').removeClass("hide");
		return;
	}

	// Switch to 'overlap' and disable 'odds ratio'
	$("#radio-value-overlaps").prop("checked", true);
	$("#radio-value-ratio-span").addClass("hide");

	// Highlight region panel
	$("#region-panel").addClass("selected");


	// Show loading animation
	$('.loader').removeClass("hide");


	getGiggleUrls();
	var dataUrl = giggleUrl + "?region=" + $('#overlaps').val();

	uploadedFileName = null;

	// get matrix data (tab delimited) and fill in heatmap
	$.ajax({
	    url: dataUrl,
	    type: "GET",
	    crossDomain: true,
	    dataType: "text",
	    success: function(data) {
	        $('.loader').addClass("hide");
	    	loadHeatmapChart(data, def);

	    },
	    error: function(error) {
	    	$('.loader').addClass("hide");
	    }
	});

}

function loadHeatmapChart(data, theDef) {
	def = theDef ? theDef : def; 
	dataForChart = data ? data : dataForChart;
	def.cells = [];

	if($("input[type='radio'].radio-value-field").is(':checked')) {
    	valueField = $("input[type='radio'].radio-value-field:checked").val();    	
	}

	dataForChart.split("\n").forEach(function(row) {
		fields = row.split("\t");
                if (fields[0] == "QueryFileError") {
                    $('#bad-file-warning').removeClass("hide");		
	            return;	
                }

		if (fields.length == 0 || fields[0] == "") {
                    
		} else {
			var rec = {};
			rec.name = fields[0].split("/")[1];
			rec.size = fields[1];
			var overlaps = fields[2];
			if (overlaps.indexOf(":") > 0) {
				rec.overlaps = overlaps.split(":")[1];
			} else {
				rec.overlaps = overlaps;
			}
			if (fields.length > 3) {
				var ratio = fields[3];
				if (ratio.indexOf(":") > 0) {
					rec.ratio = ratio.split(":")[1];
				} else {
					rec.ratio = ratio;
				}				
			}
			if (fields.length > 4) {
				var sig = fields[4];
				if (sig.indexOf(":") > 0) {
					rec.sig = sig.split(":")[1];
				} else {
					rec.sig = sig;
				}	
				rec.sig = rec.sig * 100;			
			}
			if (fields.length > 5) {
				var combo = fields[5];
				if (combo.indexOf(":") > 0) {
					rec.combo = combo.split(":")[1];
				} else {
					rec.combo = combo;
				}	
			}

			rec.row = sourceFileMap[rec.name].row;
			rec.col = sourceFileMap[rec.name].col;
			def.cells.push(rec);
		}
	});

        if ( def.cells.length > 0) {
            maxValue = d3.max(def.cells, function(d,i) {return d.overlaps});
            if (maxValue <= 2) {
                    maxValue = 3;
            }
            heatmap.colors(colorbrewer.YlGnBu[Math.min(maxValue, 9)])

            
            var selection = d3.select("#chart").datum(def);
            heatmap(selection);	
        }
}

function addCommas(nStr)
{
	nStr += '';
	x = nStr.split('.');
	x1 = x[0];
	x2 = x.length > 1 ? '.' + x[1] : '';
	var rgx = /(\d+)(\d{3})/;
	while (rgx.test(x1)) {
		x1 = x1.replace(rgx, '$1' + ',' + '$2');
	}
	return x1 + x2;
}
