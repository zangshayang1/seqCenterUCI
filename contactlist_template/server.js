var express = require('express');
var app = express();

var mongojs = require('mongojs');
var db = mongojs('contactlist', ['contactlist']);

var bodyParser = require('body-parser');

app.use(express.static(__dirname + '/public'));
// load static files such as file.html, file.css, file.js, image.file etc.
// __dirname -> current dir path

app.use(bodyParser.json());


app.get('/jobSubApp', function(req, res) {
    
    console.log("received a GET request.");
    
    // why is there a function() in find() ? 
    // that is what mongojs does so every db updates/query on server end will trigger a response to client.
    db.contactlist.find(function(err, docs) {
        // test if server receives the data from DB
        console.log(docs);
        // send back to controller
        res.json(docs); 
    });
});


// listen to the post request from the controller
app.post('/jobSubApp', function(req, res) {
    // console.log(req.body);  need another module to teach server how to parse the body
    
    // insert req.body and trigger a function to give response
    db.contactlist.insert(req.body, function(err, docs) {
        res.json(docs);
    });
});

app.listen(8888);
console.log("server running");

