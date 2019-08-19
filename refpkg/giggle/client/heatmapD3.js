function heatmapD3() {
	var dispatch = d3.dispatch("d3click");

	var margin = { top: 10, right: 100, bottom: 100, left: 5 };
	var legendMargin = 30;

	var cellSize = 20;
	var legendCellSize = 60;
	var rowLabelWidth = 420;
	var colLabelHeight = 110;

	var score    = function(d) { return +d.overlaps; };

	var rowLabels  = function(d) { return d.dimensions[0].elements};
	var colLabels  = function(d) { return d.dimensions[1].elements};
	var cellData  = function(d) { return d.cells};
	var cellRow   = function(d) { return +d.row; };
	var cellCol   = function(d) { return +d.col; };
	var cellValue = function(d) { return +d.overlaps; };


	var buckets = 9;
	var colors = ["#ffffd9","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58"]; // alternatively colorbrewer.YlGnBu[9]

	var rowNames = [];
	var colNames = [];

	var rowIndex = {};
	var colIndex = {};

	defaults = {};

	function chart(selection) {

		selection.each(function(jsonData) {

			var data = cellData(jsonData);
			rowNames = rowLabels(jsonData);
			colNames = colLabels(jsonData);

			var width  = (colNames.length*cellSize);
			var height = (rowNames.length*cellSize);

			var colorScale = d3.scale
			                   .quantize()
			                   .domain([d3.min(data, function (d) { return cellValue(d); }), buckets - 1, d3.max(data, function (d) { return cellValue(d); })])
			                   .range(colors);

			d3.select('#chart svg').remove();

			var svg = d3.select("#chart")
			            .append("svg")
			            .attr("width", width + rowLabelWidth + margin.left + legendMargin + margin.right)
			            .attr("height", height + colLabelHeight + margin.top + margin.bottom);

			var rowGroup = svg.append("g")
			                  .attr("transform", "translate(" 
			                   + (rowLabelWidth + margin.left - 10)
			                   + "," 
			                   + (+margin.top + (colLabelHeight) + (cellSize/2)) + ")");


			rowGroup.selectAll(".rowLabel")
			        .data(rowNames)
			        .enter().append("text")
			        .text(function (d) { return d; })
			        .attr("x", 0)
			        .attr("y", function (d, i) { return i * cellSize; })
			        .style("text-anchor", "end")
			        .attr("class", "rowLabel");


			var colGroup = svg.append("g")
			                  .attr("transform", "translate(" 
					                  + (margin.left + rowLabelWidth + (cellSize/3))
					                  + "," + (+margin.top + colLabelHeight) +")");

			colGroup.selectAll(".colLabel")
			      .data(colNames)
			      .enter()
			      .append("g")
			      .attr("transform", function(d,i) {
			         return "translate(" + (i * cellSize) + ",0)";
			      })
			      .append("text")
			      .attr("dx", ".8em")
			      .attr("dy", ".15em")
			      .style("text-anchor", "start")
			      .style("transform", "rotate(-45deg)")
			      .attr("class", "colLabel")
			      .text(function(d) { return d; });

			for (var i = 0; i < rowNames.length; i++) {
				var slot = rowNames[i];
				rowIndex[slot] = i;
			}
			for (var i = 0; i < colNames.length; i++) {
				var slot = colNames[i];
				colIndex[slot] = i;
			}

			var cellGroup = svg.append("g")
				               .attr("transform", "translate(" + 
				                     (rowLabelWidth + +margin.left + cellSize) + "," + (colLabelHeight + +margin.top + cellSize) + ")");

			var cells = cellGroup.selectAll(".score")
					             .data(data);

			

			cells.enter()
				 .append("rect")
				  .attr("x", function(d) { 				    
				    return (cellCol(d) - 1) * cellSize; 
				  })
				  .attr("y", function(d) { 
				    return (cellRow(d) - 1) * cellSize; 
				  })
				  .attr("rx", 0)
				  .attr("ry", 0)
				  .attr("class", "score bordered")
				  .attr("width", cellSize)
				  .attr("height", cellSize)
				  .style("fill", colors[0])
				  .on("click", function(d) {
					  	var cellObject = $().extend(d);
					  	cellObject.rowLabel = rowNames[cellRow(d)];
					  	cellObject.colLabel = colNames[cellCol(d)];
			            dispatch.d3click(cellObject);
		          });

			cells.transition()
			     .duration(1000)
			     .style("fill", function(d) { 
			        return colorScale(cellValue(d)); 
			     });

			

			cells.exit().remove();

			var interval = (colorScale.domain()[1] - colorScale.domain()[0]) / colorScale.range().length;
			var breaks = d3.range(0, colorScale.range().length).map(function(i) { return colorScale.domain()[0] + d3.round(i * interval); });

			var legendGroup = svg.append("g")
			                     .attr("transform", "translate(" + 
				                    (rowLabelWidth + +margin.left + width + legendMargin) 
				                    + "," 
				                    + (colLabelHeight + margin.top) + ")");

			legendGroup.selectAll(".legend").remove();
			var legend = legendGroup.selectAll(".legend")
			                        .data(breaks);

			legend.enter().append("g")
			  			 .attr("class", "legend");

			legend.append("rect")
					.attr("x", 0)
					.attr("y", function(d, i) { return legendCellSize * i; })
					.attr("width", cellSize / 2)
					.attr("height", legendCellSize)
					.style("fill", function(d, i) { return colors[i]; });

			legend.append("text")
					.attr("class", "legendLabel")
					.text(function(d) { return "≥ " + Math.round(d); })
					.style("text-anchor", "start")
					.attr("x", cellSize)
					.attr("y",  function(d, i) { 
					  return (legendCellSize * i) + legendCellSize/2; 
					})

			legend.exit().remove();

        }); 				
	}
	chart.cellValue = function(_) {
		if (!arguments.length) return cellValue;
			cellValue = _;
		return chart;
	};
	chart.rowLabels = function(_) {
		if (!arguments.length) return rowLabels;
			rowLabels = _;
		return chart;
	};
	chart.colLabels = function(_) {
		if (!arguments.length) return colLabels;
			colLabels = _;
		return chart;
	};
	chart.cellData = function(_) {
		if (!arguments.length) return cellData;
			cellData = _;
		return chart;
	};
	chart.cellRow = function(_) {
		if (!arguments.length) return cellRow;
			cellRow = _;
		return chart;
	};	
	chart.cellCol = function(_) {
		if (!arguments.length) return cellCol;
			cellCol = _;
		return chart;
	};	

	chart.margin = function(_) {
		if (!arguments.length) return margin;
			margin = _;
		return chart;
	};
	chart.cellSize = function(_) {
		if (!arguments.length) return cellSize;
			cellSize = _;
		return chart;
	};
	chart.legendCellSize = function(_) {
		if (!arguments.length) return legendCellSize;
			legendCellSize = _;
		return chart;
	};
	chart.rowLabelWidth = function(_) {
		if (!arguments.length) return rowLabelWidth;
			rowLabelWidth = _;
		return chart;
	};
	chart.colLabelHeight = function(_) {
		if (!arguments.length) return colLabelHeight;
			colLabelHeight = _;
		return chart;
	};
	chart.buckets = function(_) {
		if (!arguments.length) return buckets;
			buckets = _;
		return chart;
	};
	chart.colors = function(_) {
		if (!arguments.length) return colors;
			colors = _;
		return chart;
	};


	// This adds the "on" methods to our custom exports
	d3.rebind(chart, dispatch, "on");

	return chart;
}
