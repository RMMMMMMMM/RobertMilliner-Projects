using Microsoft.AspNetCore.Identity;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Threading.Tasks;



namespace Infrastructure.Models;

// Add profile data for application users by adding properties to the AppUser class
public class AppUser : IdentityUser
{

    [Required]
    [DisplayName("First Name")]
    public string? FName {  get; set; }

    [Required]
    [DisplayName("Last Name")]
    public string? LName { get; set; }

  
    [DisplayName("Street Address")]
    public string? StreetAddress { get; set; }


    [DisplayName("City")]
    public string? City { get; set; }

    [Required]
    [DisplayName("State")]
    public string? State { get; set; }

    [Required]
    [DisplayName("Postal Code")]
    public string? PostalCode { get; set;}

    [NotMapped]
    public string FullName { get { return FName + " " + LName; } }

}

