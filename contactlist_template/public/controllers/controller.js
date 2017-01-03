var myApp = angular.module('myApp',[]);

myApp.controller('appController', ['$scope', 
                                    '$http', 
                                    function($scope, 
                                             $http) {
    console.log("controller initiated and said: hello world.");
    
    // send a request to server
    // .success is what follows if the http request(get, post) is successful
    // it stores what came back from server to var response
    $http.get('/jobSubApp').success(function(response) {
    $scope.contactlist = response;
    // allow us to access these data in VIEW (index.html).

    // create an addContact function for ng-click to call
    $scope.addContact = function() {
        // post to server behavior will be triggered every time
        // "add" button is clicked
        $http.post('/jobSubApp', $scope.contact).success(function(response) {
            console.log("we successfully added:", response);
        });
    };
    });
}]);
