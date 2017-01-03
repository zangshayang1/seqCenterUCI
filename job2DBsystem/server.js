// server js

var express = require("express");
var app = express();

var mongojs = require("mongojs");
var mongodb = mongojs("seqJob", ["seqJob"]);

var bodyParser = require("body-parser");

app.use(express.static(__dirname + "/public"));
app.use(bodyParser.json());



// listen to POST request from seqApp, receive req, send res
app.post("/job2DBsystem", function(req, res) {
   mongodb.seqJob.insert(req.body, function(err, docs) {
        // insert req.body and trigger reponse function
        console.log("new entry wrote to seqJob successful!")
        res.json(docs);
    });
});

// listen to GET request from seqApp, receive req, send res 
app.get("/job2DBsystem", function(req, res) {
    mongodb.seqJob.find().sort({"_id": -1}).limit(15).toArray(function(err, docs) {
        // sort() so that it shows the top 15 latest jobs with the latest top1 coming job will be put on the top
        res.json(docs);
        console.log("GET request returned.");
    });
});

// keep listening at port 8888
app.listen(8888);
console.log("server listening...")
