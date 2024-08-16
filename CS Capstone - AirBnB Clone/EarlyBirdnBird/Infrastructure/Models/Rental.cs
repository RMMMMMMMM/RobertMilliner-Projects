using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Models
{
    public class Rental
    {
        [Key]
        public int Id { get; set; }

        //[ForeignKey("AmenityId")]
        //public Amenity? Amenity { get; set; }

        //[ForeignKey("BlockId")]
        //public Block? Block { get; set; }

        //[ForeignKey("ImageId")]
        //public Image? Image { get; set; }

        //[ForeignKey("PriceId")]
        //public Price? Price { get; set; }

        //[ForeignKey("BedId")]
        //public Bed? Bed { get; set; }

        [ForeignKey("UserName")]
        [Display(Name = "Owner")]
        public AppUser? OwnerName { get; set; }

        // -------------------------------------------- Other attributes --------------------------------

        [Required]
        public string? Name { get; set; }

        [Required]
        public string? Description { get; set; }

        [Required]
        public string? Address { get; set; }

        [Required]
        public string? City { get; set; }

        [Required]
        public string? State { get; set; }

        [Required]
        public string? Zip { get; set; }
        
        [Required]
        public string? RentalType { get; set; }


        public int? MaxOcupants { get; set; }
        public int? Bathrooms { get; set; } // number of bathrooms
        public int? Beds { get; set; } //number of bed
        public bool? Published { get; set; } //has the user publishe the listeing 





        //these are used for many to many relationships
        public ICollection<RentalAmenity>? RentalAmenities { get; set; }
        public ICollection<RentalBed>? RentalBeds { get; set; }
        public ICollection<Price>? Prices { get; set; }
        public ICollection<Image>? Images { get; set; }
        public ICollection<Review>? Reviews { get; set; }
    }
}
