var data;
var ymin_glob=1.0;
var ymax_glob=0.0;

var margin = {top: 40.5, right: 40.5, bottom: 50.5, left: 60.5},
    width = 400 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .domain([0, 1.0])
    .range([0, width]);

var y = d3.scale.linear()
    .range([height,0]);

var y2 = d3.scale.linear()
    .range([height,0]);

var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

var line = d3.svg.line()
    .x(function(d) {return x(d.time);})
    .y(function(d) {return y(d.yield);})

var line2 = d3.svg.line()
    .x(function(d) {return x(d.time);})
    .y(function(d) {return y2(d.temp);})

var svg = d3.select("body").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

svg.append("text")
    .style("font-size", "34px")
    .style("font-weight", "bold")
    .attr("transform", "translate("+ (width/12) +  ",30)")
    .text("Results");


svg.append("g")
    .attr("class", "axis axis--x")
    .attr("transform", "translate(0," + height + ")")
    .call(xAxis)
    .append("text")
    .style("font-size", "12px")
    .attr("transform", "translate("+ (width/2) + ", 30 )")
    .text("Time [s]");


var draw = function(error, pre, fit) {
    var construct  = function(data) {
        // transform data for plotting
        color.domain(d3.keys(data[0]).filter(function(key) {
                return !(key == "time" || key == "temp");
            })
        );

        // parse time seperately
        data.forEach(function(d) {
            d.time = +d.time;
            d.temp = +d.temp;
        });

        // construct yields 
        var yields = color.domain().map(function(name) {
            return { name: name,
                     values: data.map(
                        function(d) {return {time: +d.time, temp: +d.temp, yield: +d[name]};
                })
            };
        });
        return yields
    };

    // alert(yields[0].name);
    yields_fit = construct(fit);
    yields_pre = construct(pre);

    y2.domain([
        d3.min(yields_pre, function(c) {return d3.min(c.values, function(v) {return v.temp; }); }),
        d3.max(yields_pre, function(c) {return d3.max(c.values, function(v) {return v.temp; }); })
    ]);


    var species_fit = svg.selectAll(".specie.fit")
      .data(yields_fit)
      .enter().append("g")
      .attr("class", "yield");

    var species_pre = svg.selectAll(".specie.pre")
      .data(yields_pre)
      .enter().append("g")
      .attr("class", "yield");

    var ymin = d3.min(
        [ d3.min(yields_fit, function(c) {return d3.min(c.values, function(v) {return v.yield; }); }),
          d3.min(yields_pre, function(c) {return d3.min(c.values, function(v) {return v.yield; }); })]
    );

    var ymax = d3.max(
        [d3.max(yields_fit, function(c) {return d3.max(c.values, function(v) {return v.yield; }); }),
         d3.max(yields_pre, function(c) {return d3.max(c.values, function(v) {return v.yield; }); })]
    );

    y.domain([ymin,ymax]);

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var yAxisRight = d3.svg.axis()
        .scale(y2)
        .orient("right");

    svg.append("g")
        .attr("class", "axis axis--y")
        .call(yAxis)
        .append("text")
        .style("font-size", "12px")
        .attr("transform", "translate(-30," + (height/2) + ")rotate(-90)")
        .text("Yield [-]");

    svg.append("g")
        .attr("class", "axis axis--y")
        .attr("transform", "translate("+ width + ",0)")
        .call(yAxisRight)
        .append("text")
        .style("font-size", "12px")
        .attr("transform", "translate(40,"+ (height/2) + ")rotate(-90)")
        .text("Temperature [K]");

    species_pre.append("path")
            .style("stroke-dasharray", ("3, 3"))
            .attr("class", "line")
            .attr("d",function(d) {return line(d.values);});

    species_fit.append("path")
        .attr("class", "line")
        .attr("d",function(d) {return line(d.values);});
    
    species_pre.append("path")
        .attr("class", "line")
        .style("stroke", "red")
        .attr("d",function(d) {return line2(d.values);});

    species_pre.append("text")
      .datum(function(d) { return {name: d.name, value: d.values[d.values.length - 1]}; })
      .attr("transform", function(d) { return "translate(" + x(d.value.time) + "," + y(d.value.yield) + ")"; })
      .attr("x", -30)
      .attr("dy", "-0.8em")
      .text(function(d) { return d.name; });
};

var color = d3.scale.category10();

// document.getElementById("compute").onclick = function () {
    // load the data
    // d3.tsv("http://127.0.0.1:5001/static/res.tsv", function(error, data) {
    //     if (error) {
    //         return console.log("there was an error loading the data: " + error);
    //     };
    //     draw(data, "dashed");
    // });
    //
    // d3.tsv("http://127.0.0.1:5001/static/fit.tsv", function(error, data) {
    //      if (error) {
    //         return console.log("there was an error loading the data: " + error);
    //      };
    //      draw(data, "smthing");
    // });
    //
    // draw(pre,fit)
    queue()
        .defer(d3.tsv, "http://127.0.0.1:5001/static/res.tsv")
        .defer(d3.tsv, "http://127.0.0.1:5001/static/fit.tsv")
        .await(draw)

//}
